library(mpath)

# Exercise repeated parallel calls with fixed data, lambda values, and folds.
# Every call should agree with serial CV and release its worker connections.
set.seed(9001)
n <- 90L
dat <- data.frame(x = rnorm(n), z = rnorm(n))
dat$y <- rnbinom(n, mu = exp(0.5 + 0.3 * dat$x - 0.2 * dat$z), size = 2)
folds <- rep(seq_len(3), length.out = n)
grid <- c(0.4, 0.15)

foreach::registerDoSEQ()
backend_before <- c(foreach::getDoParName(), foreach::getDoParWorkers())
open_connections <- function() {
    connections <- showConnections(all = TRUE)
    rownames(connections)[connections[, "isopen"] == "opened" &
                          connections[, "class"] %in% c("sockconn", "servsockconn")]
}
before <- open_connections()
serial <- suppressWarnings(cv.glmregNB(y ~ x + z, data = dat,
    foldid = folds, nfolds = 3, lambda = grid, plot.it = FALSE,
    parallel = FALSE))
for (iteration in seq_len(if (.Platform$OS.type == "unix") 55L else 12L)) {
    parallel_fit <- suppressWarnings(cv.glmregNB(y ~ x + z, data = dat,
        foldid = folds, nfolds = 3, lambda = grid, plot.it = FALSE,
        parallel = TRUE, n.cores = 2))
    stopifnot(isTRUE(all.equal(parallel_fit$cv, serial$cv, tolerance = 1e-8)),
              isTRUE(all.equal(parallel_fit$residmat, serial$residmat,
                               tolerance = 1e-8)),
              !length(setdiff(open_connections(), before)),
              identical(c(foreach::getDoParName(), foreach::getDoParWorkers()),
                        backend_before))
}

# The full fit succeeds but one training fold has a constant response and
# fails. The error must propagate and leave no open worker sockets, and the
# next ordinary parallel CV call must still succeed.
bad <- dat
bad$y <- 1
bad$y[c(1, 4, 7)] <- c(4, 9, 2)
err <- tryCatch(suppressWarnings(cv.glmregNB(y ~ x + z, data = bad,
    foldid = folds, nfolds = 3, lambda = grid, plot.it = FALSE,
    parallel = TRUE, n.cores = 2)), error = identity)
stopifnot(inherits(err, "error"),
          grepl("constant response", conditionMessage(err)),
          !length(setdiff(open_connections(), before)),
          identical(c(foreach::getDoParName(), foreach::getDoParWorkers()),
                    backend_before))
again <- suppressWarnings(cv.glmregNB(y ~ x + z, data = dat,
    foldid = folds, nfolds = 3, lambda = grid, plot.it = FALSE,
    parallel = TRUE, n.cores = 2))
stopifnot(isTRUE(all.equal(again$cv, serial$cv, tolerance = 1e-8)),
          !length(setdiff(open_connections(), before)))

# Direct parallel fitting must likewise leave the registered backend intact.
fit <- suppressWarnings(glmregNB(y ~ x + z, data = dat, lambda = grid,
                                  parallel = TRUE, n.cores = 2))
stopifnot(all(is.finite(fit$beta)),
          identical(c(foreach::getDoParName(), foreach::getDoParWorkers()),
                    backend_before),
          !length(setdiff(open_connections(), before)))

# A fork worker that returns no result must not silently remove a CV fold.
parallel_map <- getFromNamespace(".mpath_parallel_map", "mpath")
if (.Platform$OS.type == "unix") {
    err <- tryCatch(suppressWarnings(parallel_map(1:2,
        function(i) if (i == 1L) i else NULL, n.cores = 2,
        backend = "fork")), error = identity)
    stopifnot(inherits(err, "error"),
              grepl("failed to return", conditionMessage(err)))
}

# Check the socket backend and an already registered external backend when
# local sockets are available. Some isolated check environments deny sockets.
external <- tryCatch(parallel::makeCluster(2), error = identity)
if (inherits(external, "error")) {
    message("Skipping PSOCK lifecycle tests: local sockets are unavailable")
} else {
    tryCatch({
        doParallel::registerDoParallel(external)
        external_before <- c(foreach::getDoParName(), foreach::getDoParWorkers())
        socket_before <- open_connections()
        for (iteration in seq_len(10L)) {
            value <- parallel_map(1:3, function(i) i * 2L,
                                  n.cores = 2, backend = "psock")
            stopifnot(identical(unlist(value), c(2L, 4L, 6L)),
                      !length(setdiff(open_connections(), socket_before)),
                      identical(c(foreach::getDoParName(), foreach::getDoParWorkers()),
                                external_before))
        }
        err <- tryCatch(parallel_map(1:3, function(i) {
            if (i == 2L) stop("intentional worker failure")
            i
        }, n.cores = 2, backend = "psock"), error = identity)
        stopifnot(inherits(err, "error"),
                  grepl("intentional worker failure", conditionMessage(err)),
                  !length(setdiff(open_connections(), socket_before)),
                  identical(c(foreach::getDoParName(), foreach::getDoParWorkers()),
                            external_before))
        # Workers owned by the caller still execute work after our calls.
        stopifnot(identical(unlist(parallel::parLapply(external, 1:2, identity)), 1:2))
    }, finally = {
        foreach::registerDoSEQ()
        parallel::stopCluster(external)
    })
}
