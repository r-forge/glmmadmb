\name{glmmADMB-package}
\alias{glmmADMB-package}
\alias{glmmADMB}
\docType{package}
\title{Generalized Linear Mixed Models Using AD Model Builder}
\description{
  Fits mixed-effects models to count data using binomial, Poisson, or
  negative binomial response distributions. Zero-inflated versions of
  Poisson and negative binomial distributions are supported.
}
\details{
  \emph{Main function:}
  \tabular{ll}{
    \code{\link{glmmadmb}} \tab fit a model
  }
  \emph{Methods:}
  \tabular{ll}{
    \code{\link{anova.glmmadmb}}   \tab anova\cr
    \code{\link{logLik.glmmadmb}}  \tab log likelihood\cr
    \code{\link{plot.glmmadmb}}    \tab plot\cr
    \code{\link{print.glmmadmb}}   \tab print\cr
    \code{\link{ranef.glmmadmb}}   \tab random effects\cr
    \code{\link{summary.glmmadmb}} \tab summary
  }
  \emph{Example data:}
  \tabular{ll}{
    \code{\link{epil2}} \tab seizure counts
  }
}
\note{
  If you find that the \pkg{glmmADMB} package does not do exactly what
  you want, you may want to modify the underlying ADMB \file{.tpl} files
  being called by the \code{glmmadmb} function, which can be found via
  \code{system.file("tpl","glmmadmb.tpl",package="glmmADMB")}.

}
\author{
  Hans Skaug, Dave Fournier, Anders Nielsen, Arni Magnusson, and Ben
  Bolker.
}
\seealso{\url{http://admb-project.org}.}
\keyword{models}
