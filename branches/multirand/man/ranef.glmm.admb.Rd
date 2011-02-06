\name{ranef.glmm.admb}
\alias{ranef.glmm.admb}
\title{Extract Random Effects}
\description{
  Extract random effects from fittted \code{glmm.admb} models, and
  (optionally) their standard deviations.
}
\usage{
ranef.glmm.admb(object, sd=FALSE, \dots)
}
\arguments{
  \item{object}{a fitted \code{glmm.admb} model.}
  \item{sd}{whether the standard deviations of the random effects should
    be returned, rather than point estimates.}
  \item{\dots}{further arguments (currently ignored).}
}
\details{
  The point estimates of the random effects are \sQuote{empirical Bayes}
  estimates in the sense that the hyperparameters (fixed effects,
  variance components, etc.) are fixed at their MLE's when the joint
  likelihood is maximized with respect to the random effects. However,
  uncertainty about the hyperparameters is taken into account when the
  standard deviations of the random effects are calculated.
}
\value{
  A matrix of estimated random effects (or standard deviations).
}
\seealso{
  \code{\link{glmm.admb}}.

  \code{\link{glmmADMB-package}}.
}
\keyword{models}
