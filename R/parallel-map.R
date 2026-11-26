# Run independent mpath jobs without registering a global foreach backend.
# The private backend argument also permits testing PSOCK cleanup on Unix.
.mpath_parallel_map <- function(X, FUN, n.cores,
                                backend = if (.Platform$OS.type == "unix") "fork" else "psock") {
    if (backend == "fork") {
        result <- parallel::mclapply(X, FUN, mc.cores = n.cores)
        failed <- which(vapply(result, inherits, logical(1), "try-error"))
        if (length(failed)) stop(as.character(result[[failed[1L]]]), call. = FALSE)
        if (any(vapply(result, is.null, logical(1))))
            stop("a parallel worker failed to return a result", call. = FALSE)
        result
    } else if (backend == "psock") {
        cl <- parallel::makeCluster(n.cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        # Workers started by Rscript need the caller's library paths, which
        # can differ from the environment's default R_LIBS settings.
        parallel::clusterCall(cl, function(paths) {
            .libPaths(paths)
            loadNamespace("mpath")
            NULL
        }, .libPaths())
        parallel::parLapply(cl, X, FUN)
    } else stop("unknown parallel backend", call. = FALSE)
}
