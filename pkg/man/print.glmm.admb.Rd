\name{print.glmm.admb}
\alias{print.glmm.admb}
\alias{coef.glmm.admb}
\alias{residuals.glmm.admb}
\alias{fitted.glmm.admb}
\title{Basic accessor methods for \code{glmm.admb} objects}
\description{Basic accessor methods for \code{glmm.admb} objects:
  print, retrieve coefficients, etc.}
\usage{
\method{print}{glmm.admb}(x,\dots)
\method{coef}{glmm.admb}(object,\dots)
\method{residuals}{glmm.admb}(object,\dots)
\method{fitted}{glmm.admb}(object,\dots)
}
\arguments{
  \item{x}{a fitted \code{glmm.admb} model.}
  \item{object}{a fitted \code{glmm.admb} model.}
  \item{\dots}{further arguments. If \code{sd_S_print=TRUE}, the
    covariance matrix of the random effects vector and corresponding
    standard deviations are obtained.}
}
\seealso{
  \code{\link{glmm.admb}}, \code{\link{glmmADMB-package}}.
}
\keyword{models}
