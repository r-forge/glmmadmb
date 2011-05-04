`\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.6 (4 May 2011)}{
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

\section{CHANGES IN VERSION 0.4-1 (26 Oct 2010)}{
  \itemize{
    \item Fixed executable bits in Linux.
    \item Improved backward compatibility with old versions of R.
}}
\section{CHANGES IN VERSION 0.4-0 (17 Oct 2010)}{
  \itemize{
    \item Removed unnecessary internal functions.
    \item Removed unnecessary warnings.
}}
\section{CHANGES IN VERSION 0.3-0 (28 Jul 2010)}{
  \itemize{
    \item Legacy version, now on R-Forge.
}}
