## modeled after summary.glm, print.summary.glm 
summary.glmm.admb <- function(object, ...) 
{
  ## print.glmm.admb(object, ...)
  ## calculate coef table

  ## for now, assume dispersion KNOWN
  ##  glm.nb inherits from glm, so summary.glm is used
  ##    don't allow for uncertainty from estimating theta??
  ## est.disp <- object$family %in% c("binom","poisson")
  est.disp <- FALSE
  
  coef.p <- object$b
  s.err <- object$stdbeta
  tvalue <- coef.p/s.err

  dn <- c("Estimate", "Std. Error")
  if(!est.disp) { # known dispersion
    pvalue <- 2*pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p),
                                 c(dn, "z value","Pr(>|z|)"))
  }
  ans <- c(object,
           list(coefficients=coef.table))
  class(ans) <- "summary.glmm.admb"
  ## modeled after summary.glm
  ans
}

print.summary.glmm.admb <- function(x, digits = max(3, getOption("digits") - 4),
                              symbolic.cor = x$symbolic.cor,
                              signif.stars = getOption("show.signif.stars"),
                              ...)
  {
    cat("\nCall:\n",
        paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    ## print deviance residuals?
    ## cat("Deviance Residuals: \n")
    ## if(x$df.residual > 5) {
    ##     x$deviance.resid <- quantile(x$deviance.resid,na.rm=TRUE)
    ##     names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    ## }
    ## xx <- zapsmall(x$deviance.resid, digits + 1)
    ## print.default(xx, digits=digits, na.print = "", print.gap = 2)
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                 na.print="NA", ...)
    cat("\n",x$n," total observations; ",x$q," groups (",x$group,
        ")\n",sep="")
    if (!is.null(x$S))
      cat("Random effect variance (",x$group,"): ",x$S," (std. err.: ",x$sd_S,")\n",
          sep="")
    if (!is.null(x$alpha))
      cat("Negative binomial alpha: ",x$alpha," (std. err.: ",x$sd_alpha,")\n",
        sep="")
    if (!is.null(om$pz))
      cat("Zero-inflation:",x$pz,"\n")

    cat("\nLog-likelihood:",x$loglik,"\n")
    ## offset
    ## cat("\nOffset:\n")
  }
    
