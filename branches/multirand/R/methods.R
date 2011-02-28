## other accessor methods (some trivial)
coef.glmm.admb <- function(object, ...) {
  object$b
}

ranef.glmm.admb <- function(object, ...) {
  object$U
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

VarCorr <- function(x,...) {
  UseMethod("VarCorr")
}

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

