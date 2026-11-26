library(mpath)

# The example from issue #7 previously divided by a zero penalty factor
# while computing the initial lambda sequence.
reported <- data.frame(A = c(15, 15, 11, 14, 0),
                       B = c(15, 13, 0, 15, 9),
                       C = c(13, 12, 0, 21, 11))
fit_reported <- suppressWarnings(zipath(C ~ A + B, data = reported,
    family = "poisson", penalty.factor.count = c(10, 0),
    penalty.factor.zero = c(10, 0), penalty = "enet", nlambda = 5))
stopifnot(inherits(fit_reported, "zipath"),
          all(is.finite(fit_reported$lambda.count)),
          all(is.finite(fit_reported$lambda.zero)))

# Strong unpenalized predictors affect both components. The initial lambda
# must be based on the model including them, not an intercept-only model.
set.seed(707)
n <- 600L
u <- rnorm(n)
v <- rnorm(n)
p <- 0.7 * u + rnorm(n)
q <- 0.6 * v + rnorm(n)
mu <- exp(0.3 + 0.65 * u + 0.3 * p)
prob <- plogis(-0.8 + 0.85 * v + 0.35 * q)
y <- rpois(n, mu)
y[rbinom(n, 1, prob) == 1] <- 0
dat <- data.frame(y = y, u = u, v = v, p = p, q = q)

fit <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", penalty.factor.count = c(10, 0),
    penalty.factor.zero = c(10, 0), nlambda = 8, maxit = 1000,
    maxit.em = 300, reltol = 1e-8, thresh = 1e-9))
null_fit <- pscl::zeroinfl(y ~ u | v, data = dat, dist = "poisson")
mu0 <- as.vector(predict(null_fit, type = "count"))
phi0 <- as.vector(predict(null_fit, type = "zero"))
prob0 <- phi0 + (1 - phi0) * exp(-mu0)
count_score <- ifelse(y == 0, -(1 - phi0) * mu0 * exp(-mu0) / prob0,
                      y - mu0)
zero_score <- ifelse(y == 0, phi0 * (1 - phi0) * (1 - exp(-mu0)) / prob0,
                     -phi0)
standardized <- function(x) (x - mean(x)) / sqrt(mean((x - mean(x))^2))
# The solver normalizes c(10, 0) to c(2, 0). Check the first lambda
# against independently evaluated KKT scores at the unpenalized null fit.
count_gradient <- abs(mean(standardized(p) * count_score))
zero_gradient <- abs(mean(standardized(q) * zero_score))
stopifnot(all(is.finite(fit$lambda.count)), all(fit$lambda.count >= 0),
          all(is.finite(fit$lambda.zero)), all(fit$lambda.zero >= 0),
          count_gradient <= 2 * fit$lambda.count[1] + 1e-6,
          zero_gradient <= 2 * fit$lambda.zero[1] + 1e-6,
          abs(fit$coefficients$count["p", 1]) < 1e-3,
          abs(fit$coefficients$zero["q", 1]) < 1e-3,
          abs(fit$coefficients$count["u", 1]) > 0.1,
          abs(fit$coefficients$zero["v", 1]) > 0.1,
          abs(fit$coefficients$count["u", 1] - coef(null_fit, model = "count")[["u"]]) < 0.05,
          abs(fit$coefficients$zero["v", 1] - coef(null_fit, model = "zero")[["v"]]) < 0.05)

# A small positive alpha is not ridge: the entry lambda must use that exact
# alpha rather than the 0.001 approximation reserved for alpha = 0.
small_alpha <- 0.0001
small <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", nlambda = 3,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    alpha.count = small_alpha, alpha.zero = small_alpha,
    maxit = 1000, maxit.em = 300, reltol = 1e-8, thresh = 1e-9))
stopifnot(abs(2 * small_alpha * small$lambda.count[1] - count_gradient) < 1e-5,
          abs(2 * small_alpha * small$lambda.zero[1] - zero_gradient) < 1e-5,
          abs(small$coefficients$count["p", 1]) < 1e-3,
          abs(small$coefficients$zero["q", 1]) < 1e-3)

# Explicit sequences are preserved and zero-penalty predictors remain free
# to enter even under large penalties on the other predictors.
grid <- c(5, 2, 0.5)
manual <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", lambda.count = grid, lambda.zero = grid,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    maxit = 1000, maxit.em = 300, reltol = 1e-8, thresh = 1e-9))
stopifnot(identical(unname(manual$lambda.count), grid),
          identical(unname(manual$lambda.zero), grid),
          all(is.finite(manual$coefficients$count)),
          all(is.finite(manual$coefficients$zero)),
          abs(manual$coefficients$count["u", 1]) > 0.1,
          abs(manual$coefficients$zero["v", 1]) > 0.1,
          abs(manual$coefficients$count["p", 1]) < 1e-3,
          abs(manual$coefficients$zero["q", 1]) < 1e-3)
nonactive <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", lambda.count = grid, lambda.zero = grid,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    type.path = "nonactive", maxit = 1000, maxit.em = 300,
    reltol = 1e-8, thresh = 1e-9))
stopifnot(all(is.finite(nonactive$coefficients$count)),
          all(is.finite(nonactive$coefficients$zero)),
          max(abs(nonactive$fitted.values - manual$fitted.values)) < 1e-3)

# The two components may have different unpenalized predictors.
mixed <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", nlambda = 4,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(0, 10),
    maxit.em = 100))
stopifnot(all(is.finite(mixed$lambda.count)), all(is.finite(mixed$lambda.zero)),
          all(is.finite(mixed$coefficients$count)),
          all(is.finite(mixed$coefficients$zero)))

# Zero factors in only one component must not disturb the other component.
for (which in c("count", "zero")) {
    factors_count <- if (which == "count") c(10, 0) else c(10, 10)
    factors_zero <- if (which == "zero") c(10, 0) else c(10, 10)
    separate <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
        family = "poisson", nlambda = 3,
        penalty.factor.count = factors_count,
        penalty.factor.zero = factors_zero, maxit.em = 100))
    stopifnot(all(is.finite(separate$lambda.count)),
              all(is.finite(separate$lambda.zero)),
              all(is.finite(separate$coefficients$count)),
              all(is.finite(separate$coefficients$zero)))
}

# Offsets and case weights must be used in the unpenalized starting model.
dat$off_count <- runif(n, -0.2, 0.2)
dat$off_zero <- runif(n, -0.2, 0.2)
dat$wt <- sample(1:2, n, replace = TRUE)
weighted <- suppressWarnings(zipath(
    y ~ p + u + offset(off_count) | q + v + offset(off_zero), data = dat,
    weights = wt, family = "poisson", nlambda = 5,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    maxit = 1000, maxit.em = 300, reltol = 1e-8, thresh = 1e-9))
weighted_null <- pscl::zeroinfl(y ~ u + offset(off_count) |
                                v + offset(off_zero), data = dat,
                                weights = wt, dist = "poisson")
stopifnot(all(is.finite(weighted$lambda.count)),
          all(is.finite(weighted$lambda.zero)),
          abs(weighted$coefficients$count["p", 1]) < 1e-3,
          abs(weighted$coefficients$zero["q", 1]) < 1e-3,
          abs(weighted$coefficients$count["u", 1] -
              coef(weighted_null, model = "count")[["u"]]) < 0.05,
          abs(weighted$coefficients$zero["v", 1] -
              coef(weighted_null, model = "zero")[["v"]]) < 0.05)

# The negative-binomial path also accepts mixed zero and positive factors.
set.seed(708)
dat$y_nb <- rnbinom(n, size = 2, mu = mu)
dat$y_nb[rbinom(n, 1, prob) == 1] <- 0
nb <- suppressWarnings(zipath(y_nb ~ p + u | q + v, data = dat,
    family = "negbin", init.theta = 2, theta.fixed = TRUE, nlambda = 4,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    maxit.em = 100))
stopifnot(all(is.finite(nb$lambda.count)), all(is.finite(nb$lambda.zero)),
          all(is.finite(nb$coefficients$count)),
          all(is.finite(nb$coefficients$zero)),
          abs(nb$coefficients$count["p", 1]) < 1e-3,
          abs(nb$coefficients$zero["q", 1]) < 1e-3)
# Independently maximize the nuisance model at fixed size=2 and check both
# penalized scores against their respective first lambda values.
nb_null_loss <- function(par) {
    m <- exp(par[1] + par[2] * dat$u)
    phi <- plogis(par[3] + par[4] * dat$v)
    a <- dnbinom(0, size = 2, mu = m)
    lp <- ifelse(dat$y_nb == 0, log(phi + (1 - phi) * a),
                 log1p(-phi) + dnbinom(dat$y_nb, size = 2, mu = m, log = TRUE))
    -sum(lp)
}
nb_null <- optim(c(0, 0, -1, 0), nb_null_loss, method = "BFGS",
                 control = list(reltol = 1e-11, maxit = 2000))
stopifnot(nb_null$convergence == 0)
mu_nb <- exp(nb_null$par[1] + nb_null$par[2] * dat$u)
phi_nb <- plogis(nb_null$par[3] + nb_null$par[4] * dat$v)
a_nb <- dnbinom(0, size = 2, mu = mu_nb)
d_nb <- phi_nb + (1 - phi_nb) * a_nb
count_nb <- ifelse(dat$y_nb == 0,
    -(1 - phi_nb) * a_nb * mu_nb * 2 / (2 + mu_nb) / d_nb,
    2 * (dat$y_nb - mu_nb) / (2 + mu_nb))
zero_nb <- ifelse(dat$y_nb == 0,
    phi_nb * (1 - phi_nb) * (1 - a_nb) / d_nb, -phi_nb)
stopifnot(abs(mean(standardized(p) * count_nb)) <= 2 * nb$lambda.count[1] + 1e-6,
          abs(mean(standardized(q) * zero_nb)) <= 2 * nb$lambda.zero[1] + 1e-6)

# The nuisance fit must also estimate theta when it is not fixed, taking
# case weights and offsets into account. Compare the initial lambda KKT
# thresholds with pscl's independently fitted nuisance model.
nb_weighted <- suppressWarnings(zipath(
    y_nb ~ p + u + offset(off_count) | q + v + offset(off_zero),
    data = dat, weights = wt, family = "negbin", theta.fixed = FALSE,
    nlambda = 4, penalty.factor.count = c(10, 0),
    penalty.factor.zero = c(10, 0), maxit = 1000, maxit.em = 300,
    reltol = 1e-8, thresh = 1e-9))
nb_weighted_null <- suppressWarnings(pscl::zeroinfl(
    y_nb ~ u + offset(off_count) | v + offset(off_zero), data = dat,
    weights = wt, dist = "negbin"))
weighted_theta <- nb_weighted_null$theta
mu_nb <- as.vector(predict(nb_weighted_null, type = "count"))
phi_nb <- as.vector(predict(nb_weighted_null, type = "zero"))
a_nb <- dnbinom(0, size = weighted_theta, mu = mu_nb)
d_nb <- phi_nb + (1 - phi_nb) * a_nb
count_nb <- ifelse(dat$y_nb == 0,
    -(1 - phi_nb) * a_nb * mu_nb * weighted_theta /
        (weighted_theta + mu_nb) / d_nb,
    weighted_theta * (dat$y_nb - mu_nb) / (weighted_theta + mu_nb))
zero_nb <- ifelse(dat$y_nb == 0,
    phi_nb * (1 - phi_nb) * (1 - a_nb) / d_nb, -phi_nb)
weighted_mean <- function(x) sum(dat$wt * x) / sum(dat$wt)
weighted_standardized <- function(x) {
    centered <- x - weighted_mean(x)
    centered / sqrt(weighted_mean(centered^2))
}
stopifnot(all(is.finite(nb_weighted$lambda.count)),
          all(is.finite(nb_weighted$lambda.zero)),
          abs(weighted_mean(weighted_standardized(p) * count_nb)) <=
              2 * nb_weighted$lambda.count[1] + 1e-5,
          abs(weighted_mean(weighted_standardized(q) * zero_nb)) <=
              2 * nb_weighted$lambda.zero[1] + 1e-5)

# All-zero factors use an unpenalized fit; ridge paths use the documented
# small positive alpha approximation to select an automatic initial lambda.
all_free <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", nlambda = 3,
    penalty.factor.count = c(0, 0), penalty.factor.zero = c(0, 0),
    maxit.em = 100))
ridge <- suppressWarnings(zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", nlambda = 3, alpha.count = 0, alpha.zero = 0,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    maxit.em = 100))
stopifnot(all(all_free$lambda.count == 0), all(all_free$lambda.zero == 0),
          all(is.finite(all_free$coefficients$count)),
          all(is.finite(all_free$coefficients$zero)),
          all(is.finite(ridge$lambda.count)), all(is.finite(ridge$lambda.zero)),
          all(is.finite(ridge$coefficients$count)),
          all(is.finite(ridge$coefficients$zero)))

# Invalid penalty factors are rejected in R, before any Fortran call.
expect_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
for (factors in list(c(1, -1), c(1, NA_real_), c(1, Inf), c(1))) {
    stopifnot(expect_error(zipath(y ~ p + u | q + v, data = dat,
                                  family = "poisson",
                                  penalty.factor.count = factors)),
              expect_error(zipath(y ~ p + u | q + v, data = dat,
                                  family = "poisson",
                                  penalty.factor.zero = factors)))
}

# Unsupported combinations must fail clearly rather than returning a fit
# that silently mishandles zero penalty factors.
stopifnot(expect_error(zipath(y ~ p + u | q + v, data = dat,
                              family = "poisson", link = "probit",
                              penalty.factor.count = c(10, 0))),
          expect_error(zipath(y ~ p + u | q + v, data = dat,
                              family = "poisson", intercept = FALSE,
                              penalty.factor.count = c(10, 0))))

# Cross-validation generates each fold's path and predicts with the same
# predictor columns. A supplied deterministic split keeps this reproducible.
folds <- split(seq_len(n), rep(seq_len(3), length.out = n))
cv <- suppressWarnings(cv.zipath(y ~ p + u | q + v, data = dat,
    family = "poisson", nlambda = 4, nfolds = 3, foldid = folds,
    penalty.factor.count = c(10, 0), penalty.factor.zero = c(10, 0),
    maxit.em = 100, plot.it = FALSE))
stopifnot(inherits(cv, "cv.zipath"), all(is.finite(cv$cv)),
          all(is.finite(cv$residmat)))
