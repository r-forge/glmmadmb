\name{epil2}
\alias{epil2}
\docType{data}
\title{Seizure Counts for Epileptics}
\description{
  Extended version of the 'epil' dataset of the MASS library. The three
  transformed variables \code{Visit}, \code{Base} and \code{Age} used by
  Booth et al. (2003) have been added to 'epil'.
}
\usage{data(epil2)}
\format{
  A data frame with 236 observations on the following 12 variables.
  \describe{
    \item{\code{y}}{a numeric vector}
    \item{\code{trt}}{a factor with levels \code{placebo}
      \code{progabide}}
    \item{\code{base}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{V4}}{a numeric vector}
    \item{\code{subject}}{a numeric vector}
    \item{\code{period}}{a numeric vector}
    \item{\code{lbase}}{a numeric vector}
    \item{\code{lage}}{a numeric vector}
    \item{Visit}{(rep(1:4,59)-2.5)/5}
    \item{Base}{log(base/4)}
    \item{Age}{log(age)}
  }
}
\references{
  Booth et al. (2003) Negative binomial loglinear mixed models
  STATISTICAL MODELLING 3 (3): 179-191
}
\examples{
data(epil2)
}
\keyword{datasets}
