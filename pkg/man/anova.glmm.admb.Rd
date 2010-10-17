\name{anova.glmm.admb}
\alias{anova.glmm.admb}
\title{Anova Table}
\description{Compute an anova table for a series of fitted models.}
\usage{
\method{anova}{glmm.admb}(object, \dots)
}
\arguments{
  \item{object}{a fitted \code{glmm.admb} model.}
  \item{\dots}{more fitted models, with the same random effects.}
}
\details{
  Basic anova table for testing the fixed effects (may be extended
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
  f1<-glmm.admb(y~Base*trt+Age+Visit, random=~Visit, group="subject",
                data=epil2, family="nbinom")
  f2<-glmm.admb(y~Base*trt+Age, random=~Visit, group="subject",
                data=epil2, family="nbinom")
  anova(f2,f1)
}
\keyword{models}
