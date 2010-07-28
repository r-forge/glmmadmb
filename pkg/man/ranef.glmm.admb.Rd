\name{ranef.glmm.admb}
\alias{ranef.glmm.admb}
\title{Extract Random Effects}
\description{
     Extract random effects from objects fitted with the function 'glmm.admb' 
	and (optionally) their standard deviations.
}
\usage{
ranef.glmm.admb(object, sd = FALSE)
}
\arguments{
  \item{object}{An object (a fit) returned from \code{glmm.admb}}
  \item{sd}{logical.  If 'TRUE', the standard deviations of the random effects are
                      returned, rather than point estimates.}
  \item{\dots}{Further arguments passed on}
}
\details{
The point estimates of the random effects are 'empirical Bayes' estimates in
the sense that the hyperparameters (fixed effects, variance components etc) are fixed
at their MLE's when the joint likelihood is maximized with respect to the random effects.
However, uncertainty about the hyperparameters is taken into account when the 
standard deviations of the random effects are calculated.
}

\value{
  A matrix of estimated random effects (or standard deviations).
}

\author{H. Skaug skaug@mi.uib.no, David Fournier otter@otter-rsch.com and
Anders Nielsen andersn@hawaii.edu}

\seealso{\code{\link{glmm.admb}}}

\keyword{models}% at least one, from doc/KEYWORDS
