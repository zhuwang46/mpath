library(mpath)

# A predictor present only in one observation becomes constant when that
# observation is held out.  Keep a variable predictor in each component so
# the tests also exercise the usual non-degenerate coefficient paths.
set.seed(55)
n <- 150L
x <- rnorm(n)
z <- rnorm(n)
y <- rpois(n, exp(0.3 + 0.4 * x))
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y = y, x = x, z = z,
                  rare = as.numeric(seq_len(n) == 1L), constant = 5,
                  fractional = 0.1)
train <- dat[-1L, , drop = FALSE]
grid <- c(0.2, 0.1, 0.05)

fit_path <- function(formula, data, ...) {
    zipath(formula, data = data, family = "poisson", standardize = TRUE,
           lambda.count = grid, lambda.zero = grid,
           maxit = 100, maxit.em = 50, ...)
}
finite_fit <- function(fit) {
    inherits(fit, "zipath") &&
        all(is.finite(fit$coefficients$count)) &&
        all(is.finite(fit$coefficients$zero)) &&
        all(is.finite(fit$fitted.values))
}

# Full data work both before and after subsampling; zero-variance count
# coefficients must remain in the original model-matrix order and be zero.
full <- fit_path(y ~ x + rare | 1, dat)
count <- fit_path(y ~ x + rare + constant + fractional | 1, train)
reduced <- fit_path(y ~ x | 1, train)
stopifnot(finite_fit(full), finite_fit(count), finite_fit(reduced),
          identical(rownames(count$coefficients$count),
                    c("(Intercept)", "x", "rare", "constant", "fractional")),
          all(count$coefficients$count[c("rare", "constant", "fractional"), ] == 0),
          isTRUE(all.equal(unname(count$fitted.values),
                           unname(reduced$fitted.values), tolerance = 1e-6)),
          isTRUE(all.equal(unname(predict(count, newdata = dat[1, , drop = FALSE])),
                           unname(predict(reduced, newdata = dat[1, , drop = FALSE])),
                           tolerance = 1e-6)))

# The same guarantees apply when a zero-model covariate becomes constant.
zero <- fit_path(y ~ x | z + rare + constant, train)
zero_reduced <- fit_path(y ~ x | z, train)
stopifnot(finite_fit(zero), finite_fit(zero_reduced),
          identical(rownames(zero$coefficients$zero),
                    c("(Intercept)", "z", "rare", "constant")),
          all(zero$coefficients$zero[c("rare", "constant"), ] == 0),
          isTRUE(all.equal(unname(zero$fitted.values),
                           unname(zero_reduced$fitted.values), tolerance = 1e-6)))

# If only one lambda sequence is supplied, it must remain unchanged while
# the other sequence is generated from the data.
mixed_count <- zipath(y ~ x + rare | z + rare, data = train,
                      family = "poisson", standardize = TRUE,
                      lambda.count = grid,
                      maxit = 100, maxit.em = 50)
mixed_zero <- zipath(y ~ x + rare | z + rare, data = train,
                     family = "poisson", standardize = TRUE,
                     lambda.zero = grid,
                     maxit = 100, maxit.em = 50)
stopifnot(finite_fit(mixed_count), finite_fit(mixed_zero),
          identical(unname(mixed_count$lambda.count), grid),
          identical(unname(mixed_zero$lambda.zero), grid),
          length(mixed_count$lambda.zero) == length(grid),
          length(mixed_zero$lambda.count) == length(grid),
          all(is.finite(mixed_count$lambda.zero)),
          all(is.finite(mixed_zero$lambda.count)))

# Weighted variance is zero when the only nonzero observation has zero
# weight; the feature must still be inert, even though its raw values differ.
w <- rep(1, n)
w[1] <- 0
weighted <- zipath(y ~ x + rare | z + rare, data = dat, weights = w,
                   family = "poisson", standardize = TRUE,
                   lambda.count = grid, lambda.zero = grid,
                   maxit = 100, maxit.em = 50)
stopifnot(finite_fit(weighted),
          all(weighted$coefficients$count["rare", ] == 0),
          all(weighted$coefficients$zero["rare", ] == 0))

# The first fold holds out the only rare observation. CV must retain the
# feature columns so held-out predictions have the correct dimensions.
folds <- split(seq_len(n), rep(seq_len(5), length.out = n))
cv <- cv.zipath(y ~ x + rare | z + rare, data = dat, foldid = folds,
                nfolds = 5, lambda.count = grid, lambda.zero = grid,
                standardize = TRUE, family = "poisson", plot.it = FALSE,
                maxit = 100, maxit.em = 50)
stopifnot(inherits(cv, "cv.zipath"), all(is.finite(cv$residmat)),
          all(is.finite(cv$cv)))

# A complete component consisting of constant predictors must work for both
# automatic and explicitly specified lambda paths, including zero penalty.
automatic <- zipath(y ~ constant | constant, data = train,
                    family = "poisson", standardize = TRUE,
                    nlambda = 3, maxit = 100, maxit.em = 50)
explicit <- fit_path(y ~ constant | constant, train)
unpenalized <- fit_path(y ~ constant | constant, train,
                        penalty.factor.count = 0, penalty.factor.zero = 0)
stopifnot(finite_fit(automatic), finite_fit(explicit), finite_fit(unpenalized),
          all(is.finite(automatic$lambda.count)),
          all(is.finite(automatic$lambda.zero)),
          identical(unname(explicit$lambda.count), grid),
          identical(unname(explicit$lambda.zero), grid),
          all(automatic$coefficients$count["constant", ] == 0),
          all(automatic$coefficients$zero["constant", ] == 0),
          all(explicit$coefficients$count["constant", ] == 0),
          all(explicit$coefficients$zero["constant", ] == 0),
          all(unpenalized$coefficients$count["constant", ] == 0),
          all(unpenalized$coefficients$zero["constant", ] == 0))

# A nonzero constant cannot be absorbed when the intercept is disabled.
err <- tryCatch(fit_path(y ~ x + constant | z, train, intercept = FALSE),
                error = identity)
stopifnot(inherits(err, "error"),
          grepl("constant count predictors.*intercept", conditionMessage(err)))

# Non-finite design values remain errors rather than being mistaken for a
# harmless constant column.
bad_x <- cbind(1, x = dat$x)
bad_z <- cbind(1, z = dat$z)
bad_x[1, 2] <- Inf
err <- tryCatch(zipath_fit(bad_x, bad_z, dat$y, weights = rep(1, n),
                           offsetx = rep(0, n), offsetz = rep(0, n)),
                error = identity)
stopifnot(inherits(err, "error"),
          grepl("model matrices.*finite", conditionMessage(err)))

# Constant columns are inert for each supported penalty and zero lambda in
# negative-binomial fits as well; compare with the corresponding reduced fit.
set.seed(98)
x_nb <- rnorm(180)
z_nb <- rnorm(180)
y_nb <- rnbinom(180, mu = exp(0.4 + 0.3 * x_nb), size = 2)
y_nb[rbinom(180, 1, 0.2) == 1] <- 0
nb_data <- data.frame(y = y_nb, x = x_nb, z = z_nb, constant = 5)
nb_grid <- c(0.2, 0.05, 0)
for (penalty in c("enet", "mnet", "snet")) {
    nb <- zipath(y ~ x + constant | z + constant, data = nb_data,
                 family = "negbin", init.theta = 2, theta.fixed = TRUE,
                 penalty = penalty, lambda.count = nb_grid,
                 lambda.zero = nb_grid, maxit = 1000, maxit.em = 200,
                 reltol = 1e-8, thresh = 1e-9)
    nb_reduced <- zipath(y ~ x | z, data = nb_data,
                         family = "negbin", init.theta = 2,
                         theta.fixed = TRUE, penalty = penalty,
                         lambda.count = nb_grid, lambda.zero = nb_grid,
                         maxit = 1000, maxit.em = 200,
                         reltol = 1e-8, thresh = 1e-9)
    stopifnot(finite_fit(nb), finite_fit(nb_reduced),
              all(nb$coefficients$count["constant", ] == 0),
              all(nb$coefficients$zero["constant", ] == 0),
              max(abs(nb$fitted.values - nb_reduced$fitted.values)) < 1e-4)
}
