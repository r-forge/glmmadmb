\name{admbControl}
\alias{admbControl}
\title{
  Control options for ADMB
}
\description{
  Determines options for running AD Model Builder
}
\usage{
admbControl(impSamp=0,maxfn=500,imaxfn=500,noinit=TRUE,shess=TRUE)
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
