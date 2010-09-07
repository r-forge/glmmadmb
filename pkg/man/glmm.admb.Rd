 \name{glmm.admb}
\alias{glmm.admb}

\title{Generalized Linear Mixed Models using AD Model Builder}
\description{
     Fits mixed-effects models to count data using Binomial, Poisson or
     negative binomial response distributions. Zero-inflated versions of 
     Poisson and negative binomial distributions are available.
}
\usage{ glmm.admb(fixed, random, group, data, family = "poisson", link, corStruct = "diag", 
		impSamp = 0, easyFlag = TRUE, zeroInflation = FALSE,
imaxfn = 10, save.dir= NULL,
offset)
}
\arguments{
  \item{fixed}{a two-sided linear formula object describing the
          fixed-effects part of the model, with the response on the
          left of a '~' operator and the terms, separated by '+'
          operators, on the right.}
  \item{random}{optionally, a one-sided formula object describing the
          random-effects part of the model. When 'random' is missing
	  an ordinary GLM without random effects is fitted.}
  \item{group}{a character string naming the main nesting variable.}
  \item{data}{a data frame containing the variables named in 'fixed', 'random' and 'group'.}
  \item{family}{a character string determining the response distribution: "poisson" or "nbinom".}
  \item{link}{a character string specifying the shape of the link function ("logit" or "probit") used for the "binomial" family.}
  \item{corStruct}{a character string specifying the covariance structure of the random effects vector. 
                 Two types of covariance matrices are are currently implemented: "diag" (diagonal matrix) and "full"
		 (positive definite matrix with all elements being estimated).}
  \item{impSamp}{integer.  The sample size in the importance sampling correction of the Laplace approximation 
		(impSamp=0 yields a plain Laplace approximation).}
  \item{easyFlag}{logical.  If 'TRUE', a faster but less robust optimization algorithm is employed (only "poisson" and "nbinom").}
  \item{zeroInflation}{logical.  If 'TRUE', a zero-inflated model is fitted (only "poisson" and "nbinom")}
  \item{imaxfn}{integer. Number of function evaluations used in intermediate optimization steps.}
  \item{save.dir}{If a quoted directory name is specified all the ADMB
    output files are saved there.}
  \item{offset}{optional offset for the model (see \code{\link{offset}})}
}
\details{

Currently, the "binomial" familiy only accepts Bernoully responce (0 or 1).

Parameterization of the negative binomial distribution: Var(Y) = E(Y)*(1+E(Y)/alpha).

Zero-inflation: With probability '1-pz' Y comes from a Poisson (or negative binomial) distribution,
and with probability 'pz' Y is zero (Bohning et al., 1999). Only available with "poisson" and "nbinom" response.

Parameters are estimated by maximum likelihood using the Laplace approximation to evaluate the marginal 
likelihood. When 'impSamp > 0' importance sampling is used to improve the
Laplace approximation (Skaug and Fournier, 2005).

If the message 'Proper convergence could not be reached' occurs, try to increase the parameter 'imaxfn'
and to set 'easyFlag = FALSE'.

}

\value{
  An object of class 'glmm.admb' representing the  model fit. 
  The generic function 'print' has a method to show the results of the fit.
  
  \item{b}{vector of fixed effects}
  \item{S}{covariance matrix of random effects}
  \item{alpha}{parameter in negative binomial distribution (only when 'family = "poisson"')}
  \item{pz}{Zero-inflation parameter (only when 'zeroInflation = TRUE')}
}

\references{

Bohning, D. et al (1999). The Zero-Inflated Poisson Model and the Decayed, Missing and Filled Teeth Index in Dental Epidemiology.
Journal of the Royal Statistical Society. Series A (Statistics in Society)  Vol. 162, No. 2 (1999), pp. 195-209.

Skaug and Fournier (2005). Automatic Evaluation of the Marginal Likelihood in Nonlinear Hierarchical Models. Unpublished
available from: http://bemata.imr.no/laplace.pdf

}
\author{H. Skaug skaug@mi.uib.no, David Fournier otter@otter-rsch.com and Anders Nielsen andersn@hawaii.edu}

\examples{
  data(epil2)
  glmm.admb(y~Base*trt+Age+Visit,random=~Visit,group="subject",data=epil2,family="nbinom")
}
\keyword{models}% at least one, from doc/KEYWORDS

