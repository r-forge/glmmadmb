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
    \code{\link{glmm.admb}} \tab fit a model
  }
  \emph{Methods:}
  \tabular{ll}{
    \code{\link{anova.glmm.admb}}   \tab anova\cr
    \code{\link{logLik.glmm.admb}}  \tab log likelihood\cr
    \code{\link{plot.glmm.admb}}    \tab plot\cr
    \code{\link{print.glmm.admb}}   \tab print\cr
    \code{\link{ranef.glmm.admb}}   \tab random effects\cr
    \code{\link{summary.glmm.admb}} \tab summary
  }
  \emph{Example data:}
  \tabular{ll}{
    \code{\link{epil2}} \tab seizure counts
  }
}
\note{
  If you find that the \pkg{glmmADMB} package does not do exactly what
  you want, you may want to modify the underlying ADMB \file{.tpl} files
  being called by the \code{glmm.admb} function:

  Poisson and negative binomial -
  \url{http://admb-project.googlecode.com/svn/trunk/examples/admb-re/glmmadmb/nbmm.tpl}.

  Bernoulli -
  \url{http://admb-project.googlecode.com/svn/trunk/examples/admb-re/glmmadmb/bvprobit.tpl}.
}
\author{
  Hans Skaug, Dave Fournier, Anders Nielsen, Arni Magnusson, and Ben
  Bolker.
}
\seealso{\url{http://admb-project.org}.}
\keyword{models}
