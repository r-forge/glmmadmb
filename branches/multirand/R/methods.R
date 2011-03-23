## other accessor methods (some trivial)
coef.glmm.admb <- function(object, ...) {
  object$b
}

## need to make sure that this plays nicely with
##  lme4 (S4 methods).  Not sure how.

## ranef <- function(object, ...) {
##   UseMethod("ranef")
##}

## setGeneric("ranef", function(object, ...) {
##     standardGeneric("ranef")
## })

## setMethod("ranef","glmm.admb",
##           function(object, sd=FALSE, ...) {
##    if(sd) return(object$sd_U)
##    mapply(sweep,object$U,lapply(object$S,function(z)sqrt(diag)),
##           MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
##  })

ranef.glmm.admb <- function(object, sd=FALSE, ...) {
  if(sd) return(object$sd_U)
  mapply(sweep,object$U,lapply(object$S,function(z)sqrt(diag(z))),
         MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
}

residuals.glmm.admb <- function(object, type=c("pearson", "response"), ...) {
  type <- match.arg(type)
  if (type=="response") {
    object$residuals
  } else {
    object$residuals/object$sd.est
  }
}

fitted.glmm.admb <- function(object, ...) {
  object$fitted
}

stdEr <- function(x, ...) {
  UseMethod("stdEr")
}

stdEr.glmm.admb <- function(x, ...) {
  x$stdbeta
}

vcov.glmm.admb <- function(object, ...) {
  outer(object$stdbeta,object$stdbeta)*object$corMat
}

nobs.glmm.admb <- function(object,...) {
  length(object$fitted)
}

## VarCorr <- function(x,...) {
##  UseMethod("VarCorr")
## }

VarCorr.glmm.admb <- function(x,...) {
  x$S
}

VarCorr.summary.glmm.admb <- VarCorr.glmm.admb

## FIXME:
##   needed:
##    update (for general convenience & to make drop1 work)
##    terms, extractAIC  (to make drop1 work)
##      for terms, do we want to save model frame? save_frame
##          (or save.frame or saveFrame) ?

