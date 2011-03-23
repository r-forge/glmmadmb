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

