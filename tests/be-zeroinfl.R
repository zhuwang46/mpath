library(mpath)

# A three-level factor occupies two coefficient rows but only one formula
# term.  The weak numeric term follows it in the model matrix; previously its
# coefficient index was incorrectly used as a term index.
set.seed(63)
n <- 900
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rpois(n, exp(0.4 + c(0, 0.8, 1.3)[as.integer(f)]))
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ f + x | 1, data = dat)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(inherits(selected, "zeroinfl"),
          identical(attr(selected$terms$count, "term.labels"), "f"),
          length(attr(selected$terms$zero, "term.labels")) == 0L)

# The same mismatch can occur in the zero-inflation component.
set.seed(64)
n <- 1500
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rpois(n, exp(1))
prob <- plogis(-1.7 + c(0, 1.6, 2.8)[as.integer(f)])
y[rbinom(n, 1, prob) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ 1 | f + x, data = dat)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(inherits(selected, "zeroinfl"),
          length(attr(selected$terms$count, "term.labels")) == 0L,
          identical(attr(selected$terms$zero, "term.labels"), "f"))

# An uninformative factor is removed as a whole while a strong numeric term
# remains, and increasing alpha can retain the full model.
set.seed(65)
n <- 1000
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rpois(n, exp(0.6 + 0.6 * x))
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ f + x | 1, data = dat)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
retained <- be.zeroinfl(fit, data = dat, alpha = 0.99)
stopifnot(identical(attr(selected$terms$count, "term.labels"), "x"),
          identical(attr(retained$terms$count, "term.labels"), c("f", "x")))

# One factor contrast can be weak even when the factor as a whole is highly
# significant. Selection must test the term jointly, not drop it because of
# that single contrast.
set.seed(70)
n <- 900
f <- factor(sample(letters[1:3], n, replace = TRUE))
y <- rpois(n, exp(0.5 + c(0, 0, 1.1)[as.integer(f)]))
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y, f)
fit <- pscl::zeroinfl(y ~ f | 1, data = dat)
stopifnot(summary(fit)$coefficients$count["fb", "Pr(>|z|)"] > 0.05,
          summary(fit)$coefficients$count["fc", "Pr(>|z|)"] < 0.001)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(identical(attr(selected$terms$count, "term.labels"), "f"))

# Factor elimination also works in the zero component.
set.seed(67)
n <- 1600
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rpois(n, exp(1))
y[rbinom(n, 1, plogis(-1 + 0.8 * x)) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ 1 | f + x, data = dat)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(identical(attr(selected$terms$zero, "term.labels"), "x"))

# A retained interaction keeps its lower-order terms, even when those terms
# have weak individual tests at the reference factor level.
set.seed(68)
n <- 1800
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rpois(n, exp(0.7 + c(0, 0.6, 1)[as.integer(f)] * x))
y[rbinom(n, 1, 0.2) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ f * x | 1, data = dat)
selected <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(identical(attr(selected$terms$count, "term.labels"),
                    c("f", "x", "f:x")))

# Negative-binomial models include an extra log(theta) summary row. It must
# not be treated as an explanatory term.
set.seed(66)
n <- 1500
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
y <- rnbinom(n, mu = exp(0.6 + 0.5 * x), size = 2)
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y, f, x)
fit <- pscl::zeroinfl(y ~ f + x | 1, data = dat, dist = "negbin")
selected <- be.zeroinfl(fit, data = dat, dist = "negbin", alpha = 0.05)
inferred <- be.zeroinfl(fit, data = dat, alpha = 0.05)
stopifnot(inherits(selected, "zeroinfl"),
          identical(attr(selected$terms$count, "term.labels"), "x"),
          identical(inferred$dist, "negbin"),
          identical(attr(inferred$terms$count, "term.labels"), "x"))

# Retain the fitted model's link, weights and offsets when refitting, even
# when neither the original model frame nor its model matrix was stored.
set.seed(63)
n <- 900
f <- factor(sample(letters[1:3], n, replace = TRUE))
x <- rnorm(n)
off <- runif(n, -0.2, 0.2)
off2 <- runif(n, -0.1, 0.1)
w <- sample(1:2, n, replace = TRUE)
y <- rpois(n, exp(0.4 + off + off2 + c(0, 0.8, 1.3)[as.integer(f)]))
y[rbinom(n, 1, 0.25) == 1] <- 0
dat <- data.frame(y, f, x, off, off2, w)
fit_formula_offset <- pscl::zeroinfl(y ~ f + x + offset(off) | 1,
    data = dat, weights = w, link = "probit", model = FALSE)
fit_argument_offset <- pscl::zeroinfl(y ~ f + x | 1, data = dat,
    weights = w, offset = off2, link = "probit", model = FALSE)
for (fit in list(fit_formula_offset, fit_argument_offset)) {
    stopifnot(is.null(fit$model), is.null(fit$x))
    selected <- be.zeroinfl(fit, data = dat)
    stopifnot(identical(attr(selected$terms$count, "term.labels"), "f"),
              identical(selected$link, "probit"),
              isTRUE(all.equal(selected$weights, fit$weights)),
              isTRUE(all.equal(selected$offset, fit$offset)))
}

# A model with intercepts only has nothing to eliminate.
fit <- pscl::zeroinfl(y ~ 1 | 1, data = dat)
selected <- be.zeroinfl(fit, data = dat)
stopifnot(inherits(selected, "zeroinfl"),
          length(attr(selected$terms$count, "term.labels")) == 0L,
          length(attr(selected$terms$zero, "term.labels")) == 0L)
