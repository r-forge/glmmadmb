\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}

\section{Changes in version 0.6.7.1 (2011-12-05)}{
  \itemize{
    \item SVN rev 148 (+?)
    \item added debug flag, debugging statements
    \item attempt to avoid temporary-directory collisions
    \item TPL file adds truncated Poisson and negative binomial
    \item update saved OwlModel
    \item added 'run' option to admbControl
    \item fixed predict bug (remove LHS)
    \item fix problems with calculating sd.est for models other than
  Bernoulli, Poisson, NB
  }
}

\section{Changes in version 0.6.6 (2011-11-23)}{
  \itemize{
    \item SVN rev 147
    \item add identity, cloglog links
    \item stub for gaussian family (not yet working!)
    \item use tolower() on family
    \item add (experimental) predict() function
    \item comment on easyFlag
    \item move ZI_kluge to admbControl
    \item fix name assignment bug for random effect levels
    
  }
}
  
\section{Changes in version 0.6.5 (2011-10-22)}{
  \itemize{
    \item SVN rev 142
    \item BRANCH MERGE with alpha version
    \item added mcmcControl, admbControl
    \item added shess, noinit options (poorly documented)
    \item added scapeMCMC, ggplot2 to Suggests: (the latter for vignette use)
  }
}

\section{Changes in version 0.6.4 (2011-08-30)}{
  \itemize{
    \item new binaries
  }
}

\section{Changes in version 0.6.3 (2011-08-08)}{
  \itemize{
    \item include more buildbot binaries
    \item added nbinom1 option
    \item fixed phase bugs in nbinom
  }
}

\section{Changes in version 0.6.2}{
  \itemize{
    \item include sd of zero-inflation in output and summary
    \item include buildbot binaries
    \item more information when failing to find binaries
  }
}

\section{Changes in version 0.6.1 (2011-07-08)}{
  \itemize{
    \item Allow random effects to be specified separately in 'random' (i.e. as in 'old' glmmADMB)
    \item Add fixef() as synonym for coef(), for lme4-compatibility
    \item Add 'start' argument (poorly tested)
}}

\section{Changes in version 0.6 (2011-05-04)}{
  \itemize{
    \item Changed formula to lme4-style: random effects, and offsets,
    included within the formula
    \item Allow multiple random effects
    \item Turn off minimal zero-inflation in non-ZI models (pz was
    always 0.0001)
    \item Require grouping variables to be factors
    \item New families: Binomial with N>1, Gamma, stub for Beta
    \item Allow alternative link functions (poorly tested): logit,
    probit, log
    \item new accessors: VarCorr(); ranef() [returns scaled random
    effects]; updated summary method (includes Wald tests); vcov();
    stdEr()
    \item Basic access to MCMC functionality
    \item function and class name changed from glmm.admb to glmmadmb
}}
\section{CHANGES IN VERSION 0.5-2 (2011-02-09)}{
  \itemize{
    \item SVN rev 54
    \item Added ranef(), stdEr(), and vcov() methods.
    \item Added 'type' argument to residuals().
    \item Improved summary() method.
    \item Added 'append' argument to dat_write().
    \item Added element 'corMat' to glmm.admb() output.
    \item Added attributes 'df' and 'nobs' to logLik() output.
    \item Added AIC to print() output.
    \item Added column 'logBroodSize' to Owls data.
    \item Added OwlModel example object.
    }}

\section{CHANGES IN VERSION 0.5-1 (2011-01-18)}{
  \itemize{
    \item SVN rev 32
    \item Fixed scrambled order of fitted values and predictions.
    }}

\section{CHANGES IN VERSION 0.5-0 (2011-01-17)}{
\itemize{
  \item SVN rev 30
  \item Added coef(), residuals(), and fitted() methods.
  \item Added bvprobit.tpl and nbmm.tpl ADMB source files.
  \item Added Owls example data.
  }}

\section{CHANGES IN VERSION 0.4-2 (2010-10-26)}{
  \itemize{
    \item SVN rev 23 [2010-11-05]
    \item Added MacOS executables.
    \item Added 'verbose' argument to glmm.admb().
    }}

\section{Changes in version 0.4-1 (2010-10-26)}{
  \itemize{
    \item Fixed executable bits in Linux.
    \item Improved backward compatibility with old versions of R.
}}
\section{Changes in version 0.4-0 (2010-10-17)}{
  \itemize{
    \item Removed unnecessary internal functions.
    \item Removed unnecessary warnings.
}}
\section{Changes in version 0.3-0 (2010-07-28)}{
  \itemize{
    \item Legacy version, now on R-Forge.
}}
