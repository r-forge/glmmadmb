library(glmmADMB)
## source("glmmadmb.R")

set.seed(1002)
nblock <- 10
nperblock <- 20
sd.u <- 1
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
d$eta <- with(d,0.2+0.5*x+r[f])
d$mu <- exp(d$eta)
d$y <- rpois(ntot,lambda=d$mu)

g1 <- glmmadmb(y~x+(1|f),family="poisson",data=d)
coef(g1)
VarCorr(g1)

g1M <- glmmadmb(y~x+(1|f),family="poisson",data=d,mcmc=TRUE,
                 mcmc.opts=mcmcControl(mcmc=100))

library(coda)
HPDinterval(g1M$mcmc)
summary(g1M$mcmc)
fixef(g1M)
## xyplot(as.mcmc(g1M$mcmc),layout=c(4,4),aspect="fill")
## densityplot(as.mcmc(g1M$mcmc),layout=c(4,4),aspect="fill")

### try example from simon.chamaille@cefe.cnrs.fr

nblock <- 10
indperblock <- 10
nperind <- 10
sd.u <- 1
ntot <- nblock*indperblock*
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
d$eta <- with(d,0.2+0.5*x+r[f])
d$mu <- exp(d$eta)
d$y <- rpois(ntot,lambda=d$mu)

mod.admb <- glmmadmb(formula=response~1+(1|group/ind),
                      data=somedata,
                      family="nbinom1",link="log",
                      admb.opts=admbControl(shess=FALSE,noinit=FALSE),
                      mcmc=TRUE)
