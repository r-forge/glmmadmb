\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}


\section{Changes in version 0.6.3 (2011-08-08)}{
  \itemize{
    \item include more buildbot binaries
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
    \item Allow random effects to be specified separately in 'random'
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
