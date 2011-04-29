\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}

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

\section{CHANGES IN VERSION 0.4-1 (2010-10-26)}{
\itemize{
  \item SVN rev 14
  \item Fixed executable bits in Linux.
  \item Improved backward compatibility with old versions of R.
}}



\section{CHANGES IN VERSION 0.4-0 (2010-10-17)}{
\itemize{
  \item SVN rev 10
  \item Improved help pages.
  \item Removed unnecessary internal functions.
  \item Removed unnecessary warnings.
  }}

\section{CHANGES IN VERSION 0.3-0 (2010-07-28)}{
\itemize{
  \item SVN rev 2
  \item Legacy version imported to R-Forge.
  }}
