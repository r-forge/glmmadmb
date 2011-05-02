\name{print.glmm.admb}
\alias{print.glmm.admb}
\alias{coef.glmm.admb}
\alias{residuals.glmm.admb}
\alias{fitted.glmm.admb}
\alias{vcov.glmm.admb}
\alias{VarCorr}
\alias{VarCorr.glmm.admb}
\alias{VarCorr.summary.glmm.admb}
\alias{VarCorr,glmm.admb-method}
\alias{VarCorr,summary.glmm.admb-method}
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
\method{VarCorr}{glmm.admb}(x,sigma=1,rdig=3)
\method{VarCorr}{summary.glmm.admb}(x,sigma=1,rdig=3)
}
\arguments{
  \item{x}{a fitted \code{glmm.admb} model}
  \item{object}{a fitted \code{glmm.admb} model}
  \item{type}{type of residuals to retrieve}
  \item{sigma}{for compatibility with \code{\link[nlme]{VarCorr}}:
    ignored}
  \item{rdig}{for compatibility with \code{\link[nlme]{VarCorr}}: ignored}
  \item{\dots}{further arguments. If \code{sd_S_print=TRUE}
    for \code{print.glmm.admb},
    print the covariance matrix of the random effects vector and corresponding
    standard deviations}
}
\details{
\itemize{
\item \code{coef} extracts (fixed-effect) coefficients from a
\code{glmm.admb} object, or the coefficient table (estimate, standard
error, Z-score, p-value) from a \code{summary.glmm.admb} object
\item \code{vcov} extracts the approximate variance-covariance matrix
 of the fixed-effect parameters, based on the information matrix
\item \code{stdEr} extracts the approximate standard errors of the
fixed-effect parameters, as above
\item \code{VarCorr} extracts the random effects variance-covariance matrices
}
}

\seealso{
  \code{\link{glmm.admb}}, \code{\link{anova.glmm.admb}},
  \code{\link{logLik.glmm.admb}}
}
\keyword{models}
