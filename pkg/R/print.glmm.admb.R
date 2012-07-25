print.glmmadmb <- function(x, ...)
{
  object <- x

  if(is.null(list(...)$sd_S_print))
    sd_S_print <- FALSE
  else
    sd_S_print <- list(...)$sd_S_print

  cat("\nGLMM's in R powered by AD Model Builder:\n\n")
  cat("  Family:", object$family, "\n")
  if(!is.null(object$alpha)) ## == "nbinom")
      cat("  alpha =", object$alpha, "\n")
  if(!is.null(object$link)) ##  && object$family=="nbinom")
      cat("  link =", object$link, "\n")
  if(object$zeroInflation)
      cat("  Zero inflation: p =", object$pz, "\n")
  cat("\nFixed effects:\n")
  cat("  Log-likelihood:", object$loglik, "\n")
  cat("  AIC:", AIC(object), "\n")
  cat("  Formula:", deparse(object$formula), "\n")
  print(object$b)

  ## FIXME: broken now that object no longer necessarily includes 'random'
  ##  fix by re-incorporating 'random'?
  if(!is.null(object$random))
  {
    cat("\nRandom effects:\n")
    ## cat("  Grouping factor:", object$group, "\n")
    ## cat("  Formula:", deparse(object$random), "\n")
    if(object$corStruct == "full") {
        cat("Structure: General positive-definite\n")
    } else {
        cat("Structure: Diagonal matrix\n")
    }
    print(VarCorr(object))
    if(sd_S_print)
    {
      cat("\nCovariance matrix of random effects vector (left) and corresponding standard deviations (right): \n\n")
      print(cbind(object$S,NA,NA,NA,NA,object$sd_S), na.print="")  ## FIXME: does this work for a non-4x4 matrix??
      cat("\nNote: The diagonal elements of the above left matrix are variances, NOT std's.\n")
    }
  }

  cat("\n","Number of observations: total=",x$n,", ",sep="")
  cat(paste(names(x$q),x$q,sep="=",collapse=", "),"\n")

  if(abs(object$gradloglik) >= 0.001)
    warning("Object has a large gradient component")

  invisible(NULL)
}
