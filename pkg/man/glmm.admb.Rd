\name{glmm.admb}
\alias{glmm.admb}
\title{Generalized Linear Mixed Models Using AD Model Builder}
\description{
  Fits mixed-effects models to count data using binomial, Poisson, or
  negative binomial response distributions. Zero-inflated versions of
  Poisson and negative binomial distributions are supported.
}
\usage{
glmm.admb(fixed, random, group, data, family="poisson", link,
          corStruct="diag", impSamp=0, easyFlag=TRUE,
          zeroInflation=FALSE, imaxfn=10, save.dir=NULL, offset,
          verbose=FALSE)
}
\arguments{
  \item{fixed}{a two-sided linear formula describing the fixed-effects
    part of the model. The response is on the left of a \code{~}
    operator and the terms, separated by \code{+} operators, are on the
    right.}
  \item{random}{optionally, a one-sided formula describing the
    random-effects part of the model. When \code{random} is missing, an
    ordinary GLM without random effects is fitted.}
  \item{group}{a string naming the main nesting variable.}
  \item{data}{a data frame containing the variables named in
    \code{fixed}, \code{random}, and \code{group}.}
  \item{family}{a string determining the response distribution:
    \code{"poisson"} or \code{"nbinom"}.}
  \item{link}{a string specifying the shape of the link function
    (\code{"logit"} or \code{"probit"}) used for the \code{"binomial"}
    family.}
  \item{corStruct}{a string specifying the covariance structure of the
    random effects vector. Two types of covariance matrices are
    currently implemented: \code{"diag"} (diagonal matrix) and
    \code{"full"} (positive definite matrix with all elements being
    estimated).}
  \item{impSamp}{the sample size in the importance sampling correction
    of the Laplace approximation (\code{impSamp=0} yields a plain
    Laplace approximation).}
  \item{easyFlag}{whether a faster but less robust optimization
    algorithm should be employed (only \code{"poisson"} and
    \code{"nbinom"} families).}
  \item{zeroInflation}{whether a zero-inflated model should be fitted
    (only \code{"poisson"} and \code{"nbinom"} families).}
  \item{imaxfn}{number of function evaluations used in intermediate
    optimization steps.}
  \item{save.dir}{optionally, a directory name where all ADMB output
    files are saved.}
  \item{offset}{optional offset for the model: should
be specified as the name of the offset variable (which should
be on the scale of the linear predictor, and included in the
data frame: see example below).}
  \item{verbose}{(logical) Print interim output from AD Model Builder run?}
}
\details{
  Currently, the \code{"binomial"} family only accepts Bernoulli
  response (0 or 1).

  Parameterization of the negative binomial distribution:

  \deqn{\mathrm{Var}(Y) = E(Y) \left[1 +
  \frac{E(Y)}{\alpha}\right]}{Var(Y) = E(Y) * (1 + E(Y)/alpha)}

  Zero-inflation: With probability \eqn{1-pz}, \eqn{Y} comes from a
  Poisson (or negative binomial) distribution, and with probability
  \eqn{pz}, \eqn{Y} is zero (Bohning et al. 1999). Only available with
  "poisson" and "nbinom" response.

  Parameters are estimated by maximum likelihood, using the Laplace
  approximation to evaluate the marginal likelihood. When \code{impSamp
  > 0}, importance sampling is used to improve the Laplace approximation
  (Skaug and Fournier 2006).

  If the message \code{'Proper convergence could not be reached'}
  occurs, try to increase the parameter \code{imaxfn} and to set
  \code{easyFlag=FALSE}.
}
\value{
  An object of class \code{"glmm.admb"} representing the model fit:
  \item{b}{vector of fixed effects.}
  \item{S}{covariance matrix of random effects.}
  \item{alpha}{parameter in negative binomial distribution (only when
  \code{family="poisson"}).}
  \item{pz}{Zero-inflation parameter (only when
  \code{zeroInflation=TRUE}).}
}
\references{
  Bohning, D., E. Dietz, and P. Schlattmann. 1999. The zero-inflated
  Poisson model and the decayed, missing and filled teeth index in
  dental epidemiology. Journal of the Royal Statistical Society. Series
  A (Statistics in Society) 162:195--209.

  Skaug, H.J. and D.A. Fournier. 2006. Automatic approximation of the
  marginal likelihood in non-Gaussian hierarchical models. Computational
  Statistics & Data Analysis 51:699--709.
}
\seealso{
  \code{\link{glmmADMB-package}}.
}
\examples{
  (fm <- glmm.admb(y~Base*trt+Age+Visit, random=~Visit, group="subject",
                  data=epil2, family="nbinom"))
data(Owls)
Owls$offset <- log(Owls$BroodSize) 
 om <- glmm.admb(SiblingNegotiation~FoodTreatment*SexParent,random=~1,group="Nest",
              family="nbinom",data=Owls,offset="offset")



}
\keyword{models}
