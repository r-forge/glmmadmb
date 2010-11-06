## other accessor methods (some trivial)

coef.glmm.admb <- function(object, ...) {
  object$b
}

residuals.glmm.admb <- function(object, ...) {
  object$residuals
}

fitted.glmm.admb <- function(object, ...) {
  object$fitted
}

 
