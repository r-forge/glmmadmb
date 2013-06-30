library(glmmADMB)  ## testing version
library(lme4)

load("multrand_batch.RData")


## kluge, for passing tests until I can get this sorted out
setMethod("VarCorr", signature(x="glmmadmb"), glmmADMB:::VarCorr.glmmadmb)
setMethod("VarCorr", signature(x="summary.glmmadmb"), glmmADMB:::VarCorr.glmmadmb)

sumfun <- function (x,times)
  UseMethod("sumfun")

tmpnamefun <- function(nn) {
  c(outer(c("min","mean","max"),nn,
          function(x,y) paste(y,x,sep=".")))
}
tmpsumfun <- function(x) c(min=min(x),mean=mean(x),max=max(x))

## will not work for 'old' glmmADMB (wrong ranef structure)
sumfun.glmmadmb <- function(x,times) {
  fixed <- coef(x)
  ransum <- unlist(lapply(ranef(x),
                   function(z)
                   apply(z,MARGIN=2,FUN=tmpsumfun)))
  LL <- logLik(x)
  rv <- unlist(lapply(VarCorr(x),diag))
  times <- round(times[3],2)
  mm <- c(fixed,c(rv),c(LL),ransum,times)
  rnames <- names(rv)
  names(mm) <- c(names(coef(x)),
                 paste("var(RE)",rnames,sep="."),
                 "logLik",
                 paste("U",tmpnamefun(rnames),sep="."),
                 "time")
  mm
}

sumfun.mer <- sumfun.merMod <- function(x,times) {
    ## should work for old/new lme4 ...
  fixed <- fixef(x)
  ransum <- unlist(lapply(ranef(x),
                          function(z)
                          apply(z,MARGIN=2,FUN=tmpsumfun)))
  LL <- logLik(x)
  rv <- sapply(VarCorr(x),c)
  times <- round(times[3],2)
  mm <- c(fixed,c(rv),c(LL),ransum,times)
  rnames <- names(rv)
  names(mm) <- c(names(fixef(x)),
                 paste("var(RE)",rnames,sep="."),
                 "logLik",
                 paste("U",tmpnamefun(rnames),sep="."),
                 "time")
  mm
}

sumfun2A <- function(modlist,tlist) {
  mapply(sumfun,modlist,tlist)
}

## does this work for multiple grouping variables?
sumfun(g2_GA,t2_GA)
sumfun(g2_lme4,t2_lme4)

sumfun2A(list(GA=g1_GA,lme4=g1_lme4),
        list(t1_GA,t1_lme4))

sumfun2A(list(GA=g2_GA,lme4=g2_lme4),
        list(t2_GA,t2_lme4))

