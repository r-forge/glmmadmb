"logLik.glmm.admb" <- function(object, ...)
{
  ret <- object$loglik
  class(ret) <- "logLik"

  return(ret)
}
