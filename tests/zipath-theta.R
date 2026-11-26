library(mpath)

# Simulated zero-inflated negative-binomial data exercise both the count and
# zero component; using several theta values exposes an accidental constant 1.
set.seed(1010)
n <- 400L
dat <- data.frame(x = rnorm(n), z = rnorm(n))
dat$y <- rnbinom(n, mu = exp(0.7 + 0.3 * dat$x), size = 2.5)
dat$y[rbinom(n, 1, plogis(-1 + 0.3 * dat$z)) == 1] <- 0
grid <- c(0.3, 0.1, 0.03)

fixed_fit <- function(theta, path = "nonactive", ...) {
    suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
        family = "negbin", lambda.count = grid, lambda.zero = grid,
        init.theta = theta, theta.fixed = TRUE, type.path = path,
        maxit = 100, maxit.em = 60, ...))
}
for (path in c("active", "nonactive")) for (theta in c(0.5, 2, 7)) {
    fit <- fixed_fit(theta, path)
    stopifnot(inherits(fit, "zipath"), length(fit$theta) == length(grid),
              all(is.finite(fit$theta)),
              isTRUE(all.equal(unname(fit$theta), rep(theta, length(grid)))),
              all(is.finite(fit$coefficients$count)),
              all(is.finite(fit$coefficients$zero)))
    # Independently evaluate the ZIP negative-binomial likelihood for every
    # fitted lambda to check theta is used in the model, not merely reported.
    X <- model.matrix(~ x + z, dat)
    count <- exp(X %*% fit$coefficients$count)
    zero <- plogis(X %*% fit$coefficients$zero)
    ll <- vapply(seq_along(grid), function(j) {
        positive <- dat$y > 0
        p0 <- dnbinom(0, mu = count[, j], size = theta)
        sum(log(zero[!positive, j] + (1 - zero[!positive, j]) * p0[!positive])) +
            sum(log1p(-zero[positive, j]) +
                dnbinom(dat$y[positive], mu = count[positive, j],
                        size = theta, log = TRUE))
    }, numeric(1))
    stopifnot(isTRUE(all.equal(unname(fit$loglik), ll, tolerance = 1e-7)))
}
low <- fixed_fit(0.5)
high <- fixed_fit(7)
stopifnot(max(abs(low$fitted.values - high$fitted.values)) > 1e-4)

# A fixed-size unpenalized fit should also optimize the requested NB
# likelihood, checked with a separate general-purpose optimizer.
unpenalized <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = c(0, 0), lambda.zero = c(0, 0),
    init.theta = 2, theta.fixed = TRUE, type.path = "nonactive",
    maxit = 1000, maxit.em = 500, reltol = 1e-9, thresh = 1e-9))
X <- model.matrix(~ x + z, dat)
nb_objective <- function(coef) {
    count <- exp(drop(X %*% coef[1:3]))
    zero <- plogis(drop(X %*% coef[4:6]))
    positive <- dat$y > 0
    p0 <- dnbinom(0, mu = count, size = 2)
    -mean(ifelse(positive,
                 log1p(-zero) + dnbinom(dat$y, mu = count, size = 2, log = TRUE),
                 log(zero + (1 - zero) * p0)))
}
nb_coef <- c(unpenalized$coefficients$count[, 2],
             unpenalized$coefficients$zero[, 2])
opt <- optim(nb_coef, nb_objective, method = "BFGS",
             control = list(maxit = 2000, reltol = 1e-11))
stopifnot(opt$convergence == 0L,
          abs(nb_objective(nb_coef) - opt$value) < 1e-5)

# A nonuniform penalty must be evaluated with the same factor normalization
# as the inner optimizer. Split coefficient signs for an independent L-BFGS-B
# solution of the fixed-theta elastic-net objective.
pen_count <- c(3, 1)
pen_zero <- c(2, 4)
norm_count <- pen_count * 2 / sum(pen_count)
norm_zero <- pen_zero * 2 / sum(pen_zero)
lam <- 0.01
penalized <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = c(lam, lam), lambda.zero = c(lam, lam),
    init.theta = 2, theta.fixed = TRUE, type.path = "nonactive",
    start = list(count = c(0.4, 0.15, -0.1), zero = c(-0.9, -0.1, 0.2), theta = 2),
    penalty.factor.count = pen_count, penalty.factor.zero = pen_zero,
    standardize = FALSE, maxit = 2000, maxit.em = 500,
    reltol = 1e-10, thresh = 1e-9))
bc <- penalized$coefficients$count[, 2]
bz <- penalized$coefficients$zero[, 2]
par <- c(bc[1], pmax(bc[-1], 0), pmax(-bc[-1], 0),
         bz[1], pmax(bz[-1], 0), pmax(-bz[-1], 0))
pen_objective <- function(par) {
    coef <- c(par[1], par[2:3] - par[4:5],
              par[6], par[7:8] - par[9:10])
    nb_objective(coef) + lam * (sum(norm_count * (par[2:3] + par[4:5])) +
                                sum(norm_zero * (par[7:8] + par[9:10])))
}
pen_opt <- optim(par, pen_objective, method = "L-BFGS-B",
    lower = c(-Inf, rep(0, 4), -Inf, rep(0, 4)), upper = rep(Inf, 10),
    control = list(maxit = 2000, factr = 1e4))
stopifnot(pen_opt$convergence == 0L,
          abs(pen_objective(par) - pen_opt$value) < 2e-5)

# Fixed theta is also respected with an automatically selected lambda path.
automatic <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", init.theta = 3, theta.fixed = TRUE, nlambda = 4,
    type.path = "nonactive", maxit.em = 40))
stopifnot(all(is.finite(automatic$lambda.count)),
          all(is.finite(automatic$lambda.zero)),
          isTRUE(all.equal(unname(automatic$theta), rep(3, automatic$nlambda))))

# For estimated dispersion, init.theta is a starting value, not a fixed value.
estimated <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = 5, theta.fixed = FALSE, type.path = "nonactive",
    maxit.em = 60))
stopifnot(all(is.finite(estimated$theta)), all(estimated$theta > 0),
          any(abs(estimated$theta - 5) > 1e-4))

# A caller's count and zero starts are expressed on the original predictor
# scale and must not be replaced by the default intercept-only model.
start <- list(count = c(0.4, 0.15, -0.1), zero = c(-0.9, -0.1, 0.2),
              theta = 2)
for (standardize in c(FALSE, TRUE)) {
    supplied <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
        family = "negbin", lambda.count = grid, lambda.zero = grid,
        init.theta = 2, theta.fixed = TRUE, start = start,
        standardize = standardize, type.path = "nonactive", maxit.em = 1))
    stopifnot(isTRUE(all.equal(unname(supplied$start$count), start$count,
                               tolerance = 1e-10)),
              isTRUE(all.equal(unname(supplied$start$zero), start$zero,
                               tolerance = 1e-10)),
              identical(unname(supplied$start$theta), 2))
    other <- start
    other$count <- c(1, -0.35, 0.3)
    other$zero <- c(-1.5, 0.3, -0.3)
    different <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
        family = "negbin", lambda.count = grid, lambda.zero = grid,
        init.theta = 2, theta.fixed = TRUE, start = other,
        standardize = standardize, type.path = "nonactive", maxit.em = 1))
    stopifnot(max(abs(supplied$coefficients$count - different$coefficients$count)) > 1e-6)
}

# For an estimated dispersion, a supplied init.theta takes precedence; when
# absent, start$theta provides the initial value.
estimated_from_init <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = 5, theta.fixed = FALSE, start = start,
    type.path = "nonactive", maxit.em = 1))
estimated_from_start <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    theta.fixed = FALSE, start = start, type.path = "nonactive", maxit.em = 1))
stopifnot(identical(unname(estimated_from_init$start$theta), 5),
          identical(unname(estimated_from_start$start$theta), 2))

# Without any supplied theta, initialization uses the negative-binomial null
# fit instead of silently defaulting to one.
default_start <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    theta.fixed = FALSE, type.path = "nonactive", maxit.em = 1))
null_start <- suppressWarnings(pscl::zeroinfl(y ~ 1 | 1, data = dat, dist = "negbin"))
stopifnot(isTRUE(all.equal(unname(default_start$start$theta), null_start$theta,
                           tolerance = 1e-8)))

# A supplied start remains effective on the zero-penalty path as well.
unpenalized_start <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = 2, theta.fixed = TRUE, start = start, maxit.em = 1,
    penalty.factor.count = c(1, 0), penalty.factor.zero = c(1, 0)))
stopifnot(isTRUE(all.equal(unname(unpenalized_start$start$count), start$count,
                           tolerance = 1e-10)),
          isTRUE(all.equal(unname(unpenalized_start$start$zero), start$zero,
                           tolerance = 1e-10)),
          all(unpenalized_start$theta == 2))

# Nonzero offsets and standardization must not change the predictor values at
# the supplied starting coefficients. With an unpenalized single EM step, both
# parameterizations should produce the same predictions on the original scale.
dat$off_count <- runif(n, -0.2, 0.2)
dat$off_zero <- runif(n, -0.2, 0.2)
dat$wt <- sample(1:2, n, replace = TRUE)
offset_start <- function(standardize) suppressWarnings(zipath(
    y ~ x + z + offset(off_count) | x + z + offset(off_zero), data = dat,
    weights = wt,
    family = "negbin", lambda.count = c(0, 0), lambda.zero = c(0, 0),
    init.theta = 2, theta.fixed = TRUE, start = start,
    standardize = standardize, type.path = "nonactive", maxit.em = 1,
    maxit = 5000, thresh = 1e-9))
offset_raw <- offset_start(FALSE)
offset_scaled <- offset_start(TRUE)
stopifnot(isTRUE(all.equal(unname(offset_raw$fitted.values),
                           unname(offset_scaled$fitted.values), tolerance = 1e-4)))

# Cross-validation must pass fixed theta to the full fit and each fold.
folds <- split(seq_len(n), rep(seq_len(3), length.out = n))
cv_at <- function(theta) suppressWarnings(cv.zipath(y ~ x + z | x + z,
    data = dat, family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = theta, theta.fixed = TRUE, type.path = "nonactive",
    maxit.em = 40, nfolds = 3, foldid = folds, plot.it = FALSE))
cv_low <- cv_at(0.5)
cv_high <- cv_at(7)
stopifnot(all(is.finite(cv_low$cv)), all(is.finite(cv_high$cv)),
          max(abs(cv_low$cv - cv_high$cv)) > 1e-4)

# The geometric family uses its defining size=1; Poisson fits remain valid.
geometric <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "geometric", lambda.count = grid, lambda.zero = grid,
    type.path = "nonactive", maxit.em = 40))
poisson <- suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "poisson", lambda.count = grid, lambda.zero = grid,
    type.path = "nonactive", maxit.em = 40))
stopifnot(all(is.finite(geometric$theta)), all(geometric$theta == 1),
          all(is.finite(poisson$coefficients$count)), is.null(poisson$theta))

# The fixed-dispersion parameter is scalar, finite and strictly positive.
expect_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
for (bad in list(NULL, 0, -1, NA_real_, Inf, c(1, 2), "two"))
    stopifnot(expect_error(fixed_fit(bad)))
stopifnot(expect_error(suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = 2, theta.fixed = TRUE,
    start = list(count = start$count, zero = start$zero, theta = 3)))))
# Non-logit zero-component fitting is not supported for caller-supplied
# starts by the underlying solver; reject instead of interpreting it as logit.
err <- tryCatch(suppressWarnings(zipath(y ~ x + z | x + z, data = dat,
    family = "negbin", lambda.count = grid, lambda.zero = grid,
    init.theta = 2, theta.fixed = TRUE, start = start, link = "probit")),
    error = identity)
stopifnot(inherits(err, "error"), grepl("logit", conditionMessage(err)))
