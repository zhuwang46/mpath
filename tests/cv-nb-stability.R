library(mpath)

# The inverse count link preserves large representable means, but an
# unrepresentable exponential must produce a clear error rather than Inf.
for (family in c(3L, 4L)) {
    means <- .Fortran("linkinv", n = 3L, eta = c(-1000, 0, 25),
                      family = family, mu = double(3), PACKAGE = "mpath")$mu
    stopifnot(all(is.finite(means)), means[2] == 1,
              abs(means[3] / exp(25) - 1) < 1e-12)
    err <- tryCatch(.Fortran("linkinv", n = 1L, eta = 1000,
                             family = family, mu = double(1), PACKAGE = "mpath"),
                    error = identity)
    stopifnot(inherits(err, "error"),
              grepl("overflow|non.?finite|finite mean|range|represent", conditionMessage(err),
                    ignore.case = TRUE))
}

# Large but representable means should yield the same finite deviances as the
# standard R families, including zero and positive counts.
outcome <- c(0, 1, 1e5)
mean_large <- c(1e150, 1e100, 1e75)
case_weights <- c(1, 2, 0.5)
for (family in c(3L, 4L)) {
    dev <- .Fortran("deveval", n = 3L, y = outcome, mu = mean_large,
                    theta = 2, weights = case_weights, family = family,
                    dev = 0, PACKAGE = "mpath")$dev
    ref <- if (family == 3L) sum(poisson()$dev.resids(outcome, mean_large, case_weights))
           else sum(MASS::negative.binomial(2)$dev.resids(outcome, mean_large, case_weights))
    stopifnot(is.finite(dev), abs(dev - ref) / max(1, abs(ref)) < 1e-10)
}

# A single outlying count can make the first undamped negative-binomial
# IRLS step overflow even though the constrained optimum is well behaved.
set.seed(82)
n <- 1000L
x <- cbind(outlier = as.numeric(seq_len(n) == 1L), z = rnorm(n))
y <- c(3000, rep(1, n - 1L))
lambda <- 0.001
scale_x <- sqrt(colMeans(sweep(x, 2L, colMeans(x))^2))
objective <- function(par) {
    beta <- par[2:3] - par[4:5]
    eta <- drop(par[1] + x %*% beta)
    -mean(dnbinom(y, size = 1, mu = exp(eta), log = TRUE)) +
        lambda * sum(scale_x * (par[2:3] + par[4:5]))
}
reference <- optim(rep(0, 5), objective, method = "L-BFGS-B",
                   lower = c(-Inf, rep(0, 4)), upper = rep(Inf, 5),
                   control = list(maxit = 5000, factr = 1e3))
stopifnot(reference$convergence == 0L)
fit <- glmreg_fit(x, y, family = "negbin", theta = 1, lambda = lambda,
                  alpha = 1, standardize = TRUE, maxit = 5000, thresh = 1e-9)
beta <- drop(fit$beta)
actual <- c(drop(fit$b0), pmax(beta, 0), pmax(-beta, 0))
stopifnot(all(is.finite(actual)), all(is.finite(fit$fitted.values)),
          abs(objective(actual) - reference$value) < 1e-5,
          max(abs(beta - (reference$par[2:3] - reference$par[4:5]))) < 0.05)

# A more severe outlier must not be silently reported as a successful fit
# with an arbitrarily poor objective. An explicit nonconvergence warning or
# error is also acceptable if the solver cannot obtain a stable solution.
y <- c(50000, rep(1, n - 1L))
reference <- optim(c(0, log(50000), 0, 0, 0), objective,
                   method = "L-BFGS-B", lower = c(-Inf, rep(0, 4)),
                   upper = rep(Inf, 5), control = list(maxit = 5000, factr = 1e3))
stopifnot(reference$convergence == 0L)
warnings_extreme <- character()
fit_extreme <- tryCatch(withCallingHandlers(
    glmreg_fit(x, y, family = "negbin", theta = 1, lambda = lambda,
               alpha = 1, standardize = TRUE, maxit = 5000, thresh = 1e-9),
    warning = function(w) {
        warnings_extreme <<- c(warnings_extreme, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = identity)
if (inherits(fit_extreme, "error")) {
    stopifnot(grepl("converg|diverg|numerical|overflow|step",
                    conditionMessage(fit_extreme), ignore.case = TRUE))
} else if (!all(fit_extreme$converged == 1)) {
    stopifnot(any(grepl("converg|diverg|numerical|overflow|step",
                        warnings_extreme, ignore.case = TRUE)))
} else {
    beta <- drop(fit_extreme$beta)
    actual <- c(drop(fit_extreme$b0), pmax(beta, 0), pmax(-beta, 0))
    stopifnot(all(is.finite(actual)),
              abs(objective(actual) - reference$value) < 1e-5)
}

# With an even larger count, nonuniform case weights and offsets, compare
# the full penalized objective to a separate constrained optimizer.
set.seed(829)
x <- cbind(outlier = as.numeric(seq_len(n) == 1L), z = rnorm(n))
case_weights <- sample(1:3, n, replace = TRUE)
offset <- runif(n, -0.2, 0.2)
y <- c(1e6, rep(1, n - 1L))
center <- drop(crossprod(case_weights, x)) / sum(case_weights)
scale_x <- sqrt(drop(crossprod(case_weights, sweep(x, 2, center)^2)) /
                sum(case_weights))
objective_weighted <- function(par) {
    beta <- par[2:3] - par[4:5]
    eta <- drop(par[1] + x %*% beta) + offset
    -sum(case_weights * dnbinom(y, size = 1, mu = exp(eta), log = TRUE)) /
        sum(case_weights) + lambda * sum(scale_x * (par[2:3] + par[4:5]))
}
reference <- optim(c(0, log(1e6), 0, 0, 0), objective_weighted,
                   method = "L-BFGS-B", lower = c(-Inf, rep(0, 4)),
                   control = list(maxit = 5000, factr = 1e3))
fit <- glmreg_fit(x, y, weights = case_weights, offset = offset,
                  family = "negbin", theta = 1, lambda = lambda,
                  alpha = 1, standardize = TRUE, maxit = 5000, thresh = 1e-9)
beta <- drop(fit$beta)
actual <- c(drop(fit$b0), pmax(beta, 0), pmax(-beta, 0))
stopifnot(reference$convergence == 0L, all(is.finite(fit$fitted.values)),
          abs(objective_weighted(actual) - reference$value) < 1e-5)

# Usual data should still agree with an independently fitted fixed-theta GLM.
set.seed(827)
n <- 240L
x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
y <- rnbinom(n, size = 2, mu = exp(0.3 + x %*% c(0.35, -0.25)))
reference <- glm(y ~ x, family = MASS::negative.binomial(2))
fit <- glmreg_fit(x, y, family = "negbin", theta = 2, lambda = 0,
                  standardize = FALSE, maxit = 5000, thresh = 1e-9)
stopifnot(max(abs(c(drop(fit$b0), drop(fit$beta)) - coef(reference))) < 1e-3)

# Numeric fold labels and explicit contrasts should be honored without the
# spurious contrasts/parallel-trace warnings reported in issue #8.
set.seed(828)
n <- 180L
dat <- data.frame(x = rnorm(n), group = factor(rep(c("A", "B", "C"), length.out = n)))
dat$y <- rnbinom(n, size = 2,
                 mu = exp(0.5 + 0.2 * dat$x + 0.2 * (dat$group == "C")))
foldid <- rep(rep(1:3, each = 3), length.out = n)
warn <- character()
fit_cv <- withCallingHandlers(
    cv.glmregNB(y ~ x + group, data = dat, foldid = foldid, nfolds = 3,
                lambda = c(0.4, 0.2), contrasts = list(group = "contr.sum"),
                plot.it = FALSE, parallel = FALSE, trace = FALSE),
    warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
stopifnot(all(is.finite(fit_cv$cv)),
          !any(grepl("non-list contrasts|trace argument is ignored", warn)),
          identical(fit_cv$fit$contrasts$group, "contr.sum"),
          identical(unname(sort(unlist(fit_cv$foldid))), seq_len(n)))
fit_cv2 <- suppressWarnings(cv.glmregNB(
    y ~ x + group, data = dat, foldid = foldid, nfolds = 3,
    lambda = c(0.4, 0.2), contrasts = list(group = "contr.sum"),
    plot.it = FALSE, parallel = FALSE, trace = FALSE))
stopifnot(isTRUE(all.equal(fit_cv$cv, fit_cv2$cv, tolerance = 1e-8)))

# One predictor must remain a matrix through formula processing and CV.
single <- suppressWarnings(cv.glmregNB(
    y ~ x, data = dat, foldid = foldid, nfolds = 3,
    lambda = c(0.4, 0.2), plot.it = FALSE, parallel = FALSE))
stopifnot(all(is.finite(single$cv)))

# Exercise the reported high-dispersion case through the public CV entry
# point, with nondegenerate outcomes in every training fold.
set.seed(18)
n_skew <- 300L
skew <- data.frame(x = c(8, rnorm(n_skew - 1L)), z = rnorm(n_skew),
                   y = c(3000, rnbinom(n_skew - 1L, mu = 2, size = 2)))
skew_folds <- rep(rep(1:3, each = 3), length.out = n_skew)
skew_cv <- suppressWarnings(cv.glmregNB(
    y ~ x + z, data = skew, foldid = skew_folds, nfolds = 3,
    lambda = c(0.1, 0.01), plot.it = FALSE, parallel = FALSE))
stopifnot(all(is.finite(skew_cv$cv)), all(is.finite(skew_cv$residmat)))

# If a fold's count response is constant, dispersion cannot be estimated
# from that fold. Report a meaningful error; a fixed-theta fit remains valid.
skew$y <- c(3000, rep(1, n_skew - 1L))
err <- tryCatch(suppressWarnings(cv.glmregNB(
    y ~ x + z, data = skew, foldid = skew_folds, nfolds = 3,
    lambda = c(0.1, 0.01), plot.it = FALSE, parallel = FALSE)), error = identity)
stopifnot(inherits(err, "error"), grepl("constant response", conditionMessage(err)))
skew$y <- 1
err <- tryCatch(suppressWarnings(glmregNB(y ~ x + z, data = skew,
    lambda = 0.1, nlambda = 1, parallel = FALSE)), error = identity)
stopifnot(inherits(err, "error"), grepl("constant response", conditionMessage(err)))
fixed <- suppressWarnings(glmregNB(y ~ x + z, data = skew, lambda = 0.1,
    nlambda = 1, theta.fixed = TRUE, theta0 = 2, parallel = FALSE))
stopifnot(all(is.finite(fixed$beta)), identical(unname(fixed$theta), 2))

# Formula and argument offsets are combined once. Missing model rows must not
# shift the response, weights, offsets, or fold indices during CV.
dat$off <- runif(n, -0.2, 0.2)
dat$extra <- runif(n, -0.1, 0.1)
dat$wt <- sample(1:2, n, replace = TRUE)
dat$x[4] <- NA_real_
keep <- complete.cases(dat)
cv_offset <- suppressWarnings(cv.glmregNB(
    y ~ x + group + offset(off), data = dat, weights = wt,
    offset = extra, foldid = foldid[keep], nfolds = 3,
    lambda = c(0.4, 0.2), plot.it = FALSE, parallel = FALSE))
stopifnot(all(is.finite(cv_offset$cv)),
          isTRUE(all.equal(unname(cv_offset$fit$offset),
                           unname((dat$off + dat$extra)[keep]))))

