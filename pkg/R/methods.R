## other accessor methods (some trivial)
coef.glmmadmb <- function(object, ...) {
  object$b
}

## for lme4/nlme compatibility
fixef.glmmadmb <- function(object, ...) {
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

## setMethod("ranef","glmmadmb",
##           function(object, sd=FALSE, ...) {
##    if(sd) return(object$sd_U)
##    mapply(sweep,object$U,lapply(object$S,function(z)sqrt(diag)),
##           MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
##  })

ranef.glmmadmb <- function(object, sd=FALSE, ...) {
  if(sd) return(object$sd_U)
  mapply(sweep,object$U,lapply(object$S,function(z)sqrt(diag(z))),
         MoreArgs=list(MARGIN=2,FUN="*"),SIMPLIFY=FALSE)
}

residuals.glmmadmb <- function(object, type=c("pearson", "response"), ...) {
  type <- match.arg(type)
  if (type=="response") {
    object$residuals
  } else {
    object$residuals/object$sd.est
  }
}

fitted.glmmadmb <- function(object, ...) {
  object$fitted
}

stdEr <- function(x, ...) {
  UseMethod("stdEr")
}

stdEr.glmmadmb <- function(x, ...) {
  x$stdbeta
}

vcov.glmmadmb <- function(object, ...) {
  outer(object$stdbeta,object$stdbeta)*object$corMat
}

nobs.glmmadmb <- function(object,...) {
  length(object$fitted)
}

## VarCorr <- function(x,...) {
##   UseMethod("VarCorr")
##}

## big difficulty here with nlme (S3 method, arguments x, sigma=1, rdig=3)
##  and lme4 (S4 methods, arguments x, ...)
VarCorr.glmmadmb <- function(x,sigma=1,rdig=3) {
  if (!missing(sigma) || !missing(rdig)) warning("'sigma' and 'rdig' arguments are present for compatibility only: ignored")
  x$S
}

VarCorr.summary.glmmadmb <- VarCorr.glmmadmb

## want to make this work when lme4 is loaded, too ... needs S4 method
setOldClass("glmmadmb")
setOldClass("summary.glmmadmb")
setMethod("VarCorr", signature(x="glmmadmb"), VarCorr.glmmadmb)
setMethod("VarCorr", signature(x="summary.glmmadmb"), VarCorr.glmmadmb)
## FIXME:
##   needed:
##    update (for general convenience & to make drop1 work)
##    terms, extractAIC  (to make drop1 work)
##      for terms, do we want to save model frame? save_frame
##          (or save.frame or saveFrame) ?


predict.glmmadmb <- function(object,...) {

}
