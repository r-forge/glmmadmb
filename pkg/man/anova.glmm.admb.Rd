\name{anova.glmm.admb}
\alias{anova.glmm.admb}

\title{Anova for glmm.admb objects}
\description{
  Calculate an anova table from a series of model fits 
}
\usage{
\method{anova}{glmm.admb}(object, \dots)
}
\arguments{
  \item{object}{A model fit returned from \code{glmm.admb}}
  \item{\dots}{One or more model fits with same random effects}
}
\details{
  Basic anova table for testing the fixed effects (may be extended later)
}
\value{
 An anova table
}
\author{H. Skaug skaug@mi.uib.no, David Fournier otter@otter-rsch.com and Anders Nielsen andersn@hawaii.edu}

\examples{
  data(epil2)
  f1<-glmm.admb(y~Base*trt+Age+Visit,random=~Visit,group="subject",data=epil2,family="nbinom")
  f2<-glmm.admb(y~Base*trt+Age,random=~Visit,group="subject",data=epil2,family="nbinom")
  anova(f2,f1)
}
\keyword{models}


