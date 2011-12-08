\name{admbControl}
\alias{admbControl}
\title{
  Control options for ADMB
}
\description{
  Determines options for running AD Model Builder
}
\usage{
admbControl(impSamp=0,maxfn=500,imaxfn=500,noinit=TRUE,shess=TRUE,
  run=TRUE, ZI_kluge=FALSE)
}
\arguments{
  \item{impSamp}{
    Number of importance sampling steps
  }
  \item{maxfn}{
    Maximum number of evaluation steps
  }
  \item{imaxfn}{number of function evaluations used in intermediate
    optimization steps.}
  \item{noinit}{
    See ADMB reference manual
  }
  \item{shess}{
    See ADMB reference manual
  }
  \item{ZI_kluge}{whether to retain a zero-inflation term of 0.001
    even when \code{zeroInflation=FALSE} (for backward
    compatibility/testing; deprecated)}
  \item{run}{Whether to run ADMB at all. This will write out files
    to the appropriate working directory, \emph{or} read existing
    output files in the working directory
    (\code{glmmadmb.\{std,bar,par,cor\}})
    into a \code{glmmadmb} object.}
}
\details{
  See the AD Model Builder reference manual
}
\value{
  Returns a list of options suitable for passing as the \code{admb.opts}
  argument to \code{\link{glmmadmb}}
}
\author{
  Ben Bolker
}
\examples{
admbControl(impSamp=200)
}
\keyword{misc}
