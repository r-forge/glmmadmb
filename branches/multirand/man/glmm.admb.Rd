\name{glmmadmb}
\alias{glmmadmb}
\title{Generalized Linear Mixed Models Using AD Model Builder}
\description{
  Fits mixed-effects models to count data using binomial, Poisson, or
  negative binomial response distributions. Zero-inflated versions of
  Poisson and negative binomial distributions are supported.
}
\usage{
glmmadmb(formula, data, family="poisson", link,start,
          random,
          corStruct="diag", impSamp=0, easyFlag=TRUE,
          zeroInflation=FALSE, ZI_kluge=FALSE, imaxfn=10,
          mcmc=FALSE,mcmc.opts=mcmc.control(),save.dir,
          verbose=FALSE,extra.args)
}
\arguments{
  \item{formula}{a two-sided linear formula describing the 
    part of the model. The response, which may be a numeric vector
    (for most models) or a matrix (for binomial models)
    is on the left of a \code{~} operator and the terms,
    separated by \code{+} operators, are on the
    right.  Random effects terms should be enclosed in parentheses
    \code{(x|g)}, where \code{x} represents the formula for the
    parameter(s) that vary across grouping variables and \code{g}
    represents the grouping variables; for example, \code{(1|g)} indicates
    a model where only intercepts vary across groups.
    Grouping variables must be factors.
    The model may include an offset term, which must be specified on the
    scale of the linear predictor and included in \code{data}:
    see \code{\link{offset}}.
    Alternatively, the random effects terms may be specified
  separately: see \code{random} below.}
  \item{data}{a data frame containing the variables named in
    \code{formula}.}
  \item{family}{a string determining the response distribution:
         currently implemented options are
	 \code{"poisson"}, \code{"binomial"}, \code{"nbinom"},
	 \code{"nbinom1"} (negative binomial ``type 1'': see Details)
	 \code{"gamma"}, and \code{"beta"}.}
  \item{link}{a string specifying the shape of the link function:
    currently implemented options are
    \code{"log"}, \code{"logit"}, and \code{"probit"}.}
  \item{start}{an (optional) named list specifying starting values, with
    (optional) components \code{fixed} (fixed-effect parameters: default
    0);
    \code{pz} (zero-inflation probability: default 0.02);
    \code{log_alpha} (log of dispersion parameter: default 1);
    \code{RE_sd} (standard deviations of random effects: default 0.25);
    \code{RE_cor} (off-diagonal elements of Cholesky decomposition of
    random-effect correlation matrix: default 0.0001);
    \code{u} (random effects: default 0)}
  \item{random}{a formula specifying the random effects.
    A single random effect may be specified as \code{~x|g},
    but parentheses must be used if more than one random
    effect is specified (e.g., \code{~(x1|g1)+(x2|g2)}).}
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
  \item{ZI_kluge}{whether to retain a zero-inflation term of 0.001
    even when \code{zeroInflation=FALSE} (for backward
    compatibility/testing; deprecated)}
  \item{imaxfn}{number of function evaluations used in intermediate
    optimization steps.}
  \item{mcmc}{(logical) run MCMC after fitting, based on estimated mode
    and parameter variance-covariance matrix}
  \item{mcmc.opts}{list of options for MCMC run: see \code{\link{mcmc.control}}}
  \item{save.dir}{Optionally, a directory name where all ADMB output
    files are saved. If \code{save.dir} is not specified, a temporary
    directory is created (with \code{\link{tempfile}}) and all ADMB input
    and output files are deleted after the ADMB run finishes and the results are read into R.}
  \item{verbose}{(logical) Print interim output from AD Model Builder
    run?}
  \item{extra.args}{(character) extra arguments to AD Model Builder (see
    ADMB manual)}
}
\details{
  \itemize{
    \item Parameterization of the negative binomial distribution:
    \itemize{
      \item for \code{family=="nbinom"}:
      \deqn{\mathrm{Var}(Y) = E(Y) \left[1 +
	\frac{E(Y)}{\alpha}\right]}{Var(Y) = E(Y) * (1 + E(Y)/alpha)}
      (i.e. Hardin and Hilbe NB2; \code{alpha} corresponds to the \code{size} parameter in
      \code{\link{dnbinom}})
      \item for \code{family=="nbinom1"}:
      \deqn{\mathrm{Var}(Y) = E(Y) \cdot \alpha}{Var(Y) = E(Y) * alpha}
      (i.e. Hardin and Hilbe NB1: \code{alpha} corresponds to the scale
      parameter in a quasi-Poisson model)
    }
    \item Zero-inflation: With probability \eqn{1-pz}, \eqn{Y} comes from a
    Poisson (or negative binomial) distribution, and with probability
    \eqn{pz}, \eqn{Y} is zero (Bohning et al. 1999).
    
    \item Parameters are estimated by maximum likelihood, using the Laplace
    approximation to evaluate the marginal likelihood.
    When \code{impSamp >0}, 
    importance sampling is used to improve the Laplace approximation
    (Skaug and Fournier 2006).
    
    \item If the message \code{'Proper convergence could not be reached'}
    occurs, try to increase the parameter \code{imaxfn} and to set
    \code{easyFlag=FALSE}.
  }
}
\value{
  An object of class \code{"glmmadmb"} representing the model fit:
  \item{b}{vector of fixed effects}
  \item{S}{covariance matrix of random effects}
  \item{alpha}{scale/overdispersion parameter
    (negative binomial, Gamma, beta)}
  \item{pz}{Zero-inflation parameter (only when
    \code{zeroInflation=TRUE})}
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
  epil2$subject <- factor(epil2$subject)
  (fm <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                  data=epil2, family="nbinom"))
\dontrun{
data(Owls)
Owls$offset <- log(Owls$BroodSize) 
 om <- glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+
           (1|Nest)+offset(offset),
## FIXME: would like to be able to use offset(log(BroodSize))
          zeroInflation=TRUE,family="nbinom",data=Owls)
}
}
\keyword{models}
