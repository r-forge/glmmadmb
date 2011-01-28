\name{print.glmm.admb}
\alias{print.glmm.admb}
\alias{coef.glmm.admb}
\alias{residuals.glmm.admb}
\alias{fitted.glmm.admb}
\alias{vcov.glmm.admb}
\alias{stdEr}
\alias{stdEr.glmm.admb}
\title{Basic accessor methods for \code{glmm.admb} objects}
\description{Basic accessor methods for \code{glmm.admb} objects:
  print, retrieve coefficients, etc.}
\usage{
stdEr(x,\ldots)
\method{print}{glmm.admb}(x,\dots)
\method{coef}{glmm.admb}(object,\dots)
\method{residuals}{glmm.admb}(object,type=c("pearson", "response"),\dots)
\method{fitted}{glmm.admb}(object,\dots)
\method{vcov}{glmm.admb}(object,\dots)
\method{stdEr}{glmm.admb}(x,\dots)
}
\arguments{
  \item{x}{a fitted \code{glmm.admb} model}
  \item{object}{a fitted \code{glmm.admb} model}
  \item{type}{type of residuals to retrieve}
  \item{\dots}{further arguments. If \code{sd_S_print=TRUE}
    for \code{print.glmm.admb},
    print the covariance matrix of the random effects vector and corresponding
    standard deviations}
}
\seealso{
  \code{\link{glmm.admb}}, \code{\link{anova.glmm.admb}},
  \code{\link{logLik.glmm.admb}}
}
\keyword{models}
