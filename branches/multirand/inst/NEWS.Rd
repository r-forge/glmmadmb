\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}

\section{Changes in version 0.6 (2011-05-04)}{
  \itemize{
    \item Changed formula to lme4-style: random effects, and offsets,
    included within the formula
    \item Allow multiple random effects
    \item Turn off minimal zero-inflation in non-ZI models
    \item Require grouping variables to be factors
    \item Allow binomial responses with N>1
    \item Allow alternative link functions (poorly tested)
    \item new accessors: VarCorr(); ranef() [returns scaled random
  effects]; updated summary method (includes Wald tests); vcov(); stdEr()
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
