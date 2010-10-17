ranef.glmm.admb <- function(object, sd=FALSE, ...)
{
  out <- if(sd) object$sd_U else object$U

  return(out)
}
