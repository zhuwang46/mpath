cv.glmregNB <- function(formula, data, weights, offset=NULL, 
lambda=NULL, nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, trace=FALSE, parallel=FALSE, contrasts=NULL,
...){
    call <- match.call()
 mf <- Call <- match.call()
                                     m <- match(c("formula", "data", "subset", "weights", "na.action",
                                                "etastart", "mustart", "offset"), names(mf), 0)
                                            mf <- mf[c(1, m)]
                                            mf$drop.unused.levels <- TRUE
                                            mf[[1L]] <- as.name("model.frame")
                                            mf <- eval.parent(mf)
  # Fit each fold using the same rows and order as the model frame, after
  # subset/NA handling. Otherwise fold indices refer to different cases.
  if (is.data.frame(data)) {
    data_rows <- match(rownames(mf), rownames(data))
    if (anyNA(data_rows)) stop("cannot align the model frame with data")
    data <- data[data_rows, , drop=FALSE]
  }
  dots <- list(...)
  # The model frame has already applied subset; do not apply it again to
  # the full or fold-specific data.
  dots$subset <- NULL
                                            Terms <- attr(mf, "terms")
                                            Y <- model.response(mf, "numeric")
  ## null model support
                                            X <- if (!is.empty.model(Terms)) model.matrix(Terms, mf, contrasts) else matrix(,NROW(Y),0)
  nobs <- n <- length(Y)
  nvars <- m <- dim(X) - 1
 weights <- model.weights(mf)
 offset <- as.vector(model.offset(mf))
 explicit_offset <- mf[["(offset)"]]
 if(!length(weights)) weights <- rep(1, nrow(mf))
  if(any(weights < 0)) stop("negative weights not allowed")

  if(missing(foldid) && nfolds < 3)
    stop("smallest nfolds should be 3\n")
  K <- nfolds
  glmregNB.obj <- do.call("glmregNB", c(list(formula, data, weights,
                                              offset=explicit_offset, lambda=lambda,
                                              contrasts=contrasts,
                                              parallel=FALSE), dots))
  lambda <- glmregNB.obj$lambda
  nlambda <- length(lambda)
  theta0 <- glmregNB.obj$theta
  if(is.null(offset)) offset <- rep(0, nobs)
  if(missing(foldid)) all.folds <- cv.folds(n, K)
  else if (is.list(foldid)) {
    all.folds <- foldid
    K <- length(all.folds)
    indices <- unlist(all.folds, use.names=FALSE)
    if (K < 3L || any(lengths(all.folds) == 0L) ||
        !is.numeric(indices) || anyNA(indices) ||
        length(indices) != n || !identical(sort(as.double(indices)), as.double(seq_len(n))))
      stop("foldid must contain at least three nonempty folds that partition the observations")
  } else {
    if (length(foldid) != n || !is.numeric(foldid) || anyNA(foldid) ||
        any(!is.finite(foldid)) || any(foldid != floor(foldid)) ||
        any(foldid < 1L))
      stop("foldid must be a vector of positive integer fold labels, one per observation")
    ids <- sort(unique(foldid))
    if (length(ids) < 3L)
      stop("foldid must identify at least three folds")
    all.folds <- lapply(ids, function(id) which(foldid == id))
    K <- length(all.folds)
  }

  foldenv <- list2env(list(formula=formula, data=data, X=X, Y=Y,
                           weights=weights, offset=offset, explicit_offset=explicit_offset, lambda=lambda,
                           nlambda=nlambda, theta0=theta0,
                           all.folds=all.folds, dots=dots, trace=trace,
                           contrasts=contrasts),
                      parent=asNamespace("mpath"))
  foldwork <- function(i) {
    omit <- all.folds[[i]]
    offsetnow <- if (is.null(explicit_offset)) NULL else explicit_offset[-omit]
    newoffset <- if (is.null(offset)) NULL else offset[omit]
    fitcv <- do.call("glmregNB", c(list(formula, data[-omit,], weights[-omit],
                                         offset=offsetnow, nlambda=nlambda,
                                         lambda=lambda, theta.est=FALSE,
                                         theta0=theta0, parallel=FALSE,
                                         trace=trace, contrasts=contrasts), dots))
    fitcv$terms <- NULL
    logLik(fitcv, newx=X[omit,-1, drop=FALSE], Y[omit],
           weights=weights[omit], offset=newoffset)
  }
  environment(foldwork) <- foldenv

  if(parallel){
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    registerDoParallel(cl)
    i <- 1
    fold.results <- foreach(i=seq_len(K), .packages="mpath") %dopar% foldwork(i)
  } else {
    fold.results <- vector("list", K)
    for(i in seq_len(K)) {
      if(trace) cat("\n CV Fold", i, "\n\n")
      fold.results[[i]] <- foldwork(i)
    }
  }
  nresid <- max(c(nlambda, lengths(fold.results)))
  residmat <- matrix(NA_real_, nresid, K)
  for(i in seq_len(K)) {
    res <- fold.results[[i]]
    if(!is.null(res)) residmat[seq_along(res), i] <- res
  }
   cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=glmregNB.obj, residmat=residmat, lambda = lambda, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which])
    class(obj) <- "cv.glmreg"
    if(plot.it) plot(obj,se=se)
    obj
  }

