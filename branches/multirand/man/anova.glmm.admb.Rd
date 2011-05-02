\name{anova.glmm.admb}
\alias{anova.glmm.admb}
\title{Model Comparison Table}
\description{Compute a model comparison (likelihood ratio test) table for a series of fitted models.}
\usage{
\method{anova}{glmm.admb}(object, \dots)
}
\arguments{
  \item{object}{a fitted \code{glmm.admb} model.}
  \item{\dots}{more fitted models, with the same random effects.}
}
\details{
  Basic table for testing the fixed effects, using
  likelihood ratio tests (may be extended
  later).
}
\value{
  An object of class \code{"anova"}, inheriting from class
  \code{"data.frame"}.
}
\seealso{
  \code{\link{glmm.admb}}.

  \code{\link{glmmADMB-package}}.
}
\examples{
  data(epil2)
  epil2$subject <- factor(epil2$subject)
  f1<-glmm.admb(y~Base*trt+Age+Visit+(Visit|subject),
                data=epil2, family="nbinom")
  f2<-glmm.admb(y~Base*trt+Age+(Visit|subject),
                data=epil2, family="nbinom")
  anova(f2,f1)
}
\keyword{models}
