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

  if(!is.null(object$random))
  {
    cat("\nRandom effects:\n")
    ## cat("  Grouping factor:", object$group, "\n")
    ## cat("  Formula:", deparse(object$random), "\n")
    if(object$corStruct == "full")
    {
      cat("Structure: General positive-definite\n")
      tmp <- cov2cor(object$S)
      tmp[upper.tri(tmp,diag=TRUE)] <- NA
      tmp <- cbind(sqrt(diag(object$S)), tmp)
      dimnames(tmp) <- list(rownames(object$S), c("StdDev","Corr",rep("",nrow(tmp)-1)))
      print(tmp, na.print="")
    }
    else
    {
      cat("Structure: Diagonal matrix\n")
      tmp <- sqrt(diag(object$S))
      names(tmp) <- colnames(object$S)
      print(tmp)
    }
    if(sd_S_print)
    {
      cat("\nCovariance matrix of random effects vector (left) and corresponding standard deviations (right): \n\n")
      print(cbind(object$S,NA,NA,NA,NA,object$sd_S), na.print="")  ## FIXME: does this work for a non-4x4 matrix??
      cat("\nNote: The diagonal elements of the above left matrix are variances, NOT std's.\n")
    }
  }

  cat("\nNumber of Observations:", object$n, "\n")
  cat("Number of Groups:", object$q, "\n")
  if(abs(object$gradloglik) >= 0.001)
    warning("Object has a large gradient component")

  invisible(NULL)
}
