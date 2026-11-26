### Backward elimination for zero-inflated regression.

# Group model-matrix columns by formula term. A factor or interaction may
# generate several columns, so its significance is a joint Wald test.
.zeroinfl_term_pvalues <- function(fit, component, data) {
    terms <- fit$terms[[component]]
    labels <- attr(terms, "term.labels")
    result <- setNames(rep(NA_real_, length(labels)), labels)
    if (!length(labels)) return(result)

    if (!is.null(fit$x) || !is.null(fit$model)) {
        mm <- model.matrix(fit, model = component)
    } else {
        # Offsets are not model-matrix columns. Rebuild without them because
        # the stored terms environment may contain an unevaluated offset arg.
        matrix_formula <- reformulate(labels, intercept = attr(terms, "intercept"),
                                      env = environment(formula(fit)))
        variables <- attr(terms(matrix_formula), "variables")[-1L]
        variables <- vapply(variables, function(x) paste(deparse(x), collapse = ""), "")
        xlev <- fit$levels[intersect(names(fit$levels), variables)]
        mf <- model.frame(matrix_formula, data = data, xlev = xlev,
                          na.action = na.pass)
        mm <- model.matrix(matrix_formula, data = mf,
                           contrasts.arg = fit$contrasts[[component]])
    }
    assignment <- attr(mm, "assign")
    beta <- fit$coefficients[[component]]
    covariance <- fit$vcov
    for (i in seq_along(labels)) {
        columns <- colnames(mm)[assignment == i]
        idx <- match(columns, names(beta))
        covnames <- paste0(component, "_", columns)
        if (!length(idx) || anyNA(idx) ||
            !all(covnames %in% rownames(covariance))) next
        b <- beta[idx]
        V <- covariance[covnames, covnames, drop = FALSE]
        if (any(!is.finite(b)) || any(!is.finite(V))) next
        statistic <- tryCatch(as.numeric(crossprod(b, solve(V, b))),
                              error = function(e) NA_real_)
        if (is.finite(statistic) && statistic >= 0)
            result[i] <- pchisq(statistic, df = length(b), lower.tail = FALSE)
    }
    result
}

be.zeroinfl <- function(object, data, dist=c("poisson", "negbin", "geometric"), alpha=0.05, trace=FALSE){
    if(!inherits(object, "zeroinfl")) stop("object must be zeroinfl\n")
    if (missing(dist)) dist <- object$dist else dist <- match.arg(dist)
    fit <- object
    history <- list()
    if (trace) {
        cat("Initial model\n")
        print(summary(fit))
    }

    repeat {
        count <- .zeroinfl_term_pvalues(fit, "count", data)
        zero <- .zeroinfl_term_pvalues(fit, "zero", data)
        # Respect marginality: do not remove a term needed by an interaction.
        count[!names(count) %in% drop.scope(fit$terms$count)] <- NA_real_
        zero[!names(zero) %in% drop.scope(fit$terms$zero)] <- NA_real_
        candidates <- c(count, zero)
        candidates[!is.finite(candidates)] <- NA_real_
        if (!any(!is.na(candidates))) break
        # Preserve the historical preference for the zero component on ties.
        best <- which(candidates == max(candidates, na.rm = TRUE))
        best <- best[length(best)]
        if (candidates[best] <= alpha) break
        component <- if (best <= length(count)) "count" else "zero"
        term_index <- if (component == "count") best else best - length(count)
        term <- names(candidates)[best]
        if (trace) {
            cat("\nstep", length(history) + 1L, "\n")
            cat("drop variable in", component, "component:", term, "\n")
        }

        count_terms <- fit$terms$count
        zero_terms <- fit$terms$zero
        if (component == "count")
            count_terms <- drop.terms(count_terms, term_index, keep.response = TRUE)
        else
            zero_terms <- drop.terms(zero_terms, term_index, keep.response = FALSE)
        count_formula <- formula(count_terms)
        zero_formula <- formula(zero_terms)
        new_formula <- as.formula(call("~", count_formula[[2L]],
                                       call("|", count_formula[[3L]], zero_formula[[2L]])),
                                  env = environment(formula(fit)))
        # Retain the original call's link, offsets, weights, and controls.
        refit_call <- fit$call
        refit_call$formula <- new_formula
        refit_call$data <- quote(data)
        refit_call$dist <- dist
        refit_env <- new.env(parent = environment(new_formula))
        refit_env$data <- data
        refit_env$zeroinfl <- zeroinfl
        next_fit <- tryCatch(eval(refit_call, envir = refit_env),
                             error = function(e) {
                                 warning("Could not refit after removing ", term,
                                         " from the ", component, " component: ",
                                         conditionMessage(e), call. = FALSE)
                                 NULL
                             })
        if (is.null(next_fit)) break
        fit <- next_fit
        history[[length(history) + 1L]] <- c(loglik = as.numeric(logLik(fit)),
                                             BIC = AIC(fit, k = log(fit$n)),
                                             AIC = AIC(fit))
        if (trace) {
            print(summary(fit))
            cat("\nloglik of zero-inflated model", logLik(fit), "\n")
            cat("\nBIC of zero-inflated model", AIC(fit, k = log(fit$n)))
            cat("\nAIC of zero-inflated model", AIC(fit), "\n")
        }
    }
    if (trace && length(history)) print(do.call(rbind, history))
    fit
}
