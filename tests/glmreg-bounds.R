library(mpath)

# Compare bounded ridge fits with an independent box-constrained optimizer.
# Correlated predictors and unequal scales make post-hoc clipping inadequate.
set.seed(423)
n <- 180L
u <- rnorm(n)
x <- cbind(0.4 * u, 2 * (0.7 * u + rnorm(n)),
           1.3 * (0.5 * u + rnorm(n)))
colnames(x) <- c("first", "second", "third")
w <- sample(1:3, n, replace = TRUE)
eta <- drop(0.2 + x %*% c(0.8, -0.35, 0.55))
responses <- list(gaussian = eta + rnorm(n, sd = 0.5),
                  binomial = rbinom(n, 1, plogis(eta)),
                  poisson = rpois(n, exp(eta)),
                  negbin = rnbinom(n, mu = exp(eta), size = 2))
off <- runif(n, -0.4, 0.4)
lower <- c(-0.25, 0, -0.3)
upper <- c(0, 0.12, 0.25)
pf <- c(1, 0.5, 2)
normalized_pf <- pf * length(pf) / sum(pf)

for (family in names(responses)) for (standardize in c(FALSE, TRUE)) {
    y <- responses[[family]]
    lambda <- if (family == "gaussian") 0.15 else if (standardize) 0.08 else 10
    center <- colSums(x * w) / sum(w)
    scales <- if (standardize)
        sqrt(colSums(sweep(x, 2, center)^2 * w) / sum(w)) else rep(1, ncol(x))
    divisor <- if (family == "gaussian" || standardize) sum(w) else 1
    objective <- function(par) {
        mu <- drop(par[1] + x %*% par[-1]) + off
        loss <- switch(family,
                       gaussian = (y - mu)^2 / 2,
                       binomial = log1p(exp(mu)) - y * mu,
                       poisson = exp(mu) - y * mu,
                       negbin = (y + 2) * log(2 + exp(mu)) - y * mu)
        sum(w * loss) / divisor +
            lambda / 2 * sum(normalized_pf * (scales * par[-1])^2)
    }
    gradient <- function(par) {
        mu <- drop(par[1] + x %*% par[-1]) + off
        residual <- switch(family, gaussian = mu - y,
                           binomial = plogis(mu) - y, poisson = exp(mu) - y,
                           negbin = 2 * (exp(mu) - y) / (2 + exp(mu)))
        c(sum(w * residual) / divisor,
          drop(crossprod(x, w * residual)) / divisor +
              lambda * normalized_pf * scales^2 * par[-1])
    }
    initial <- c(0, rep(0, ncol(x)))
    reference <- optim(initial, objective, gradient, method = "L-BFGS-B",
                       lower = c(-Inf, lower), upper = c(Inf, upper),
                       control = list(factr = 1e4, pgtol = 1e-9, maxit = 3000))
    stopifnot(reference$convergence == 0L)
    fit <- glmreg_fit(x, y, weights = w, family = family, lambda = lambda,
                      alpha = 0, penalty.factor = pf, standardize = standardize,
                      lower.limits = lower, upper.limits = upper,
                      offset = off, theta = 2, maxit = 5000, thresh = 1e-10)
    actual <- c(drop(fit$b0), drop(fit$beta))
    stopifnot(all(actual[-1] >= lower - 1e-8),
              all(actual[-1] <= upper + 1e-8),
              abs(objective(actual) - reference$value) <
                  2e-5 * max(1, abs(reference$value)),
              max(abs(actual - reference$par)) < 2e-3)

}

# Clamping an unconstrained answer cannot refit correlated coefficients.
set.seed(91)
a <- rnorm(200)
b <- 0.8 * a + rnorm(200, sd = 0.35)
xx <- cbind(a, b)
yy <- 2 * a + 0.1 * b + rnorm(200, sd = 0.2)
loss <- function(par) mean((yy - drop(par[1] + xx %*% par[-1]))^2) / 2 +
    0.005 * sum(par[-1]^2)
free <- optim(c(0, 0, 0), loss, method = "BFGS",
              control = list(reltol = 1e-12, maxit = 1000))$par
clamped <- c(free[1], 0, free[3])
clamped[1] <- optimize(function(b0) loss(c(b0, clamped[-1])),
                       interval = c(-10, 10))$minimum
constrained <- glmreg_fit(xx, yy, lambda = 0.01, alpha = 0,
                          standardize = FALSE, upper.limits = c(0, Inf),
                          maxit = 5000, thresh = 1e-10)
stopifnot(loss(clamped) > loss(c(drop(constrained$b0),
                                drop(constrained$beta))) + 1e-5)

# For elastic net, split each coefficient into nonnegative positive and
# negative parts, making the nonsmooth L1 objective independently optimizable
# with box constraints.
for (standardize in c(FALSE, TRUE)) {
    y <- responses$gaussian
    lambda <- 0.035
    alpha <- 0.6
    center <- colSums(x * w) / sum(w)
    scales <- if (standardize)
        sqrt(colSums(sweep(x, 2, center)^2 * w) / sum(w)) else rep(1, ncol(x))
    objective <- function(par) {
        b <- par[2:4] - par[5:7]
        residual <- drop(par[1] + x %*% b) - y
        sum(w * residual^2) / (2 * sum(w)) + lambda * sum(normalized_pf *
            (alpha * scales * (par[2:4] + par[5:7]) +
             (1 - alpha) / 2 * (scales * b)^2))
    }
    gradient <- function(par) {
        b <- par[2:4] - par[5:7]
        residual <- drop(par[1] + x %*% b) - y
        gb <- drop(crossprod(x, w * residual)) / sum(w) +
            lambda * normalized_pf * (1 - alpha) * scales^2 * b
        l1 <- lambda * normalized_pf * alpha * scales
        c(sum(w * residual) / sum(w), gb + l1, -gb + l1)
    }
    reference <- optim(rep(0, 7), objective, gradient, method = "L-BFGS-B",
                       lower = c(-Inf, rep(0, 6)),
                       upper = c(Inf, upper, -lower),
                       control = list(factr = 1e4, pgtol = 1e-9, maxit = 3000))
    stopifnot(reference$convergence == 0L)
    fit <- glmreg_fit(x, y, weights = w, lambda = lambda, alpha = alpha,
                      penalty.factor = pf, standardize = standardize,
                      lower.limits = lower, upper.limits = upper,
                      maxit = 5000, thresh = 1e-10)
    b <- drop(fit$beta)
    actual <- c(drop(fit$b0), pmax(b, 0), pmax(-b, 0))
    stopifnot(abs(objective(actual) - reference$value) < 2e-5,
              max(abs(b - (reference$par[2:4] - reference$par[5:7]))) < 2e-3)
}

# Infinite defaults leave the existing unconstrained path unchanged.
y <- responses$gaussian
grid <- c(0.3, 0.1, 0.02)
plain <- glmreg_fit(x, y, weights = w, lambda = grid, alpha = 0.5)
explicit <- glmreg_fit(x, y, weights = w, lambda = grid, alpha = 0.5,
                       lower.limits = -Inf, upper.limits = Inf)
stopifnot(isTRUE(all.equal(plain$beta, explicit$beta, tolerance = 1e-10)),
          isTRUE(all.equal(plain$b0, explicit$b0, tolerance = 1e-10)))

# Scalars, zero-fixed coefficients, formula interface and cross-validation.
fixed <- glmreg_fit(x, y, lambda = grid, lower.limits = c(0, -Inf, -Inf),
                    upper.limits = c(0, Inf, Inf))
stopifnot(all(fixed$beta[1, ] == 0))
pos <- glmreg(x, y, lambda = grid, lower.limits = 0)
stopifnot(all(pos$beta >= -1e-10))
dat <- data.frame(y, x)
form <- glmreg(y ~ first + second + third, data = dat, lambda = grid,
               upper.limits = 0)
stopifnot(all(form$beta <= 1e-10))
set.seed(2)
cv <- cv.glmreg(x, y, lambda = grid, nfolds = 3, plot.it = FALSE,
                lower.limits = lower, upper.limits = upper)
stopifnot(all(cv$fit$beta >= lower - 1e-10),
          all(cv$fit$beta <= upper + 1e-10), all(is.finite(cv$cv)))
set.seed(3)
cv_one <- cv.glmreg(y ~ first, data = dat, nlambda = 4, nfolds = 3,
                    lower.limits = 0, plot.it = FALSE)
stopifnot(length(cv_one$lambda) > 1L, all(is.finite(cv_one$lambda)),
          all(cv_one$fit$beta >= -1e-10), all(is.finite(cv_one$cv)))
# Automatic lambda generation accounts for unpenalized variables.
mixed <- glmreg_fit(x[, 1:2], y, nlambda = 5,
                    penalty.factor = c(0, 1), lower.limits = c(-Inf, 0))
stopifnot(length(mixed$lambda) == 5, all(is.finite(mixed$lambda)),
          all(mixed$beta[2, ] >= -1e-10))
mixed_explicit <- glmreg_fit(x[, 1:2], y, lambda = 0.1,
                             penalty.factor = c(0, 1), lower.limits = 0)
stopifnot(all(mixed_explicit$beta >= -1e-10))

expect_error <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(expect_error(glmreg_fit(x, y, lambda = 0.1, lower.limits = 0.1)),
          expect_error(glmreg_fit(x, y, lambda = 0.1, upper.limits = -0.1)),
          expect_error(glmreg_fit(x, y, lambda = 0.1, lower.limits = c(0, -1))),
          expect_error(glmreg_fit(x, y, lambda = 0.1, upper.limits = NA_real_)),
          expect_error(glmreg_fit(x, y, lambda = 0.1, lower.limits = 0,
                                  intercept = FALSE)),
          expect_error(glmreg_fit(x, y, lambda = 0.1, lower.limits = 0,
                                  penalty = "mnet")),
          expect_error(glmreg_fit(x, y, lambda = 0.1, upper.limits = 0,
                                  penalty = "snet")),
          expect_error(glmreg_fit(x[, 1:2], y, nlambda = 5,
                                  penalty.factor = c(0, 1),
                                  lower.limits = c(0, -Inf))))
