cv.zipath <- function(x, ...) UseMethod("cv.zipath", x)
 
cv.zipath.default <- function(x, ...) {
     if (extends(class(x), "Matrix"))
         return(zipath.matrix(x=x, ...))
     stop("no method for objects of class ", sQuote(class(x)),
          " implemented")
 }
 
cv.zipath.formula <- function(formula, data, weights, offset=NULL,                 contrast=NULL, ...){
     ## call and formula
     cl <- match.call()
     if(!attr(terms(formula), "intercept"))
         stop("non-intercept model is not implemented")
     if(missing(data)) data <- environment(formula)
     mf <- match.call(expand.dots = FALSE)
     m <- match(c("formula", "data", "na.action", "weights", "offset"), names(mf), 0)
     mf <- mf[c(1, m)]
     mf$drop.unused.levels <- TRUE
 
     ## extended formula processing
     if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
     {
         ff <- formula
         formula[[3]][1] <- call("+")
         mf$formula <- formula
         ffc <- . ~ .
         ffz <- ~ .
         ffc[[2]] <- ff[[2]]
         ffc[[3]] <- ff[[3]][[2]]
         ffz[[3]] <- ff[[3]][[3]]
         ffz[[2]] <- NULL
     } else {
         ffz <- ffc <- ff <- formula
         ffz[[2]] <- NULL
     }
     if(inherits(try(terms(ffz), silent = TRUE), "try-error")) {
         ffz <- eval(parse(text = sprintf( paste("%s -", deparse(ffc[[2]])), deparse(ffz) )))
     }
 
     ## call model.frame()
     mf[[1]] <- as.name("model.frame")
     mf <- eval(mf, parent.frame())
 
     ## extract terms, model matrices, response
     mt <- attr(mf, "terms")
     mtX <- terms(ffc, data = data)
     X <- model.matrix(mtX, mf)
     Xold <- X
 ### This means X possibly has intercept (1st) column depending on the specifications 
mtZ <- terms(ffz, data = data)
     mtZ <- terms(update(mtZ, ~ .), data = data)
     Z <- model.matrix(mtZ, mf)
     Zold <- Z
     Y <- model.response(mf, "numeric")
     ## weights and offset
     n <- length(Y)
     weights <- model.weights(mf)
     offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
     if(is.null(offsetx)) offsetx <- 0
   if(length(offsetx) == 1) offsetx <- rep.int(offsetx, n)
   offsetx <- as.vector(offsetx)
   offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
   if(is.null(offsetz)) offsetz <- 0
   if(length(offsetz) == 1) offsetz <- rep.int(offsetz, n)
   offsetz <- as.vector(offsetz)
     RET <- cv.zipath_fit(X, Z, Y, data, weights=weights, offsetx=offsetx, offsetz=offsetz, ...)
     RET$call <- match.call()
     RET

}
cv.zipath.matrix <- function(X, Z, Y, data, weights, offsetx=NULL, offsetz=NULL, ...){
     RET <- cv.zipath_fit(X, Z, Y, data, weights, offsetx=offsetx, offsetz=offsetz, ...)
     RET$call <- match.call()
     return(RET)
 }

cv.zipath_fit <- function(X, Z, Y, data, weights, offsetx, offsetz, family=c("poisson", "negbin", "geometric"), nlambda=100, lambda.count=NULL, 
		      lambda.zero=NULL, nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, 
              trace=FALSE, parallel=FALSE, ...){
    call <- match.call()
    family <- match.arg(family)
    if(missing(foldid) && nfolds < 3)
        stop("smallest nfolds should be 3\n")
    n <- dim(data)[1]
    #    nm <- dim(data)
#    nobs <- n <- nm[1]
#    nvars <- m <- nm[2]
#    mf <- Call <- match.call()
#    m <- match(c("formula", "data", "subset", "weights", "na.action",
#                 "offset"), names(mf), 0)
#    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
#    mf[[1L]] <- as.name("model.frame")
#    Y <- data[,all.vars(terms(formula))[[1]]]
    ## null model support
#    weights <- model.weights(mf)
#    if(!length(weights)) weights <- rep(1, length(Y))
#    if(any(weights < 0)) stop("negative weights not allowed")

    K <- nfolds
    tim <- proc.time()
    zipath.obj <- zipath_fit(X, Z, Y, weights, offsetx, offsetz, nlambda=nlambda, lambda.count=lambda.count, lambda.zero=lambda.zero, ...)
    #zipath.obj <- do.call("zipath", list(formula, data, weights, nlambda=nlambda, lambda.count=lambda.count, lambda.zero=lambda.zero, ...))
    cat("zipath_fit time within cv.zipath\n")
    print(proc.time()-tim)
    #lambda.count <- zipath.obj$lambda.count
    #lambda.zero <- zipath.obj$lambda.zero
    #nlambda <- zipath.obj$nlambda
    if(missing(foldid))
        all.folds <- cv.folds(n, K)
    else {
        all.folds <- foldid
        K <- nfolds <- length(foldid)
    }
    bic <- matrix(NA, nlambda, K)
    if(parallel){
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind, .packages='mpath') %dopar% {
        omit <- all.folds[[i]]
    #	fitcv <- do.call("zipath", list(formula, data[-omit,], weights[-omit], lambda.count=lambda.count, lambda.zero=lambda.zero, nlambda=nlambda, ...))
        fitcv <- zipath_fit(X[-omit,,drop=FALSE], Z[-omit,,drop=FALSE], Y[-omit], weights[-omit], offsetx[-omit], offsetz[-omit], nlambda=nlambda, lambda.count=lambda.count, lambda.zero=lambda.zero, ...)
        logLik(fitcv, X=X[omit,, drop=FALSE], Z=Z[omit,, drop=FALSE], Y[omit], weights=weights[omit], offsetx=offsetx[omit], offsetxz=offsetz[omit])
    }
    stopImplicitCluster()
    }
    else{
     residmat <- matrix(NA, nlambda, K)
     for(i in seq(K)) {
         cat("\n CV Fold", i, "\n\n")
       if(trace)
         cat("\n CV Fold", i, "\n\n")
       omit <- all.folds[[i]]
       fitcv <- zipath_fit(X[-omit,,drop=FALSE], Z[-omit,,drop=FALSE], Y[-omit], weights[-omit], offsetx[-omit], offsetz[-omit], nlambda=nlambda, lambda.count=lambda.count, lambda.zero=lambda.zero, ...)
       #fitcv <- do.call("zipath", list(formula, data[-omit,], weights[-omit], lambda.count=lambda.count, lambda.zero=lambda.zero, nlambda=nlambda, ...))
       residmat[, i] <- logLik(fitcv, X=X[omit,, drop=FALSE], Z=Z[omit,, drop=FALSE], Y[omit], weights=weights[omit], offsetx=offsetx[omit], offsetz=offsetz[omit])
     }
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=zipath.obj, residmat=residmat, bic=bic, cv = cv, cv.error = cv.error, foldid=all.folds, nlambda=nlambda, lambda.which= lambda.which, lambda.optim = list(count=lambda.count[lambda.which], zero=lambda.zero[lambda.which]))
    class(obj) <- c("cv.zipath")
    if(plot.it) plot(obj,se=se)
    obj
}

"plot.cv.zipath" <-
    function(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...){
        nlambda <- x$nlambda
        cv <- x$cv
        cv.error <- x$cv.error
        if(is.null(ylab))
            ylab <- "log-likelihood"
        plot(1:nlambda, cv, type = "b", xlab = "Index of lambda pair", ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
if(se)
    error.bars(1:nlambda, cv + cv.error, cv - cv.error,
               width = width, col=col)

invisible()
}

predict.cv.zipath <- function(object, newdata, ...)
predict(object$fit, newdata, which=object$lambda.which, ...)

coef.cv.zipath <- function(object, which=object$lambda.which, model = c("full", "count", "zero"), ...) {
    model <- match.arg(model)
    rval <- object$fit$coefficients
    rval <- switch(model,
                   "full" = list(count=rval$count[, which], zero=rval$zero[, which]),
                   "count" = rval$count[,which],
                   "zero" = rval$zero[,which])
    rval
}
