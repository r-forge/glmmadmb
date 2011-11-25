## testing identity, cloglog links, gaussian family
library(glmmADMB)

set.seed(1002)
nblock <- 10
nperblock <- 50
sd.u <- 1
beta <- c(0.2,0.5)
ntot <- nblock*nperblock
d <- data.frame(x=runif(ntot),f=factor(rep(LETTERS[1:nblock],each=nperblock)))
r <- rnorm(nblock,mean=0,sd=sd.u)
gshape <- 1.5
d$offset <- rgamma(ntot,1,1)
d <- within(d,
            {
              eta0 <- beta[1]+beta[2]*x+offset
              eta <- eta0+r[f]
            })


## cloglog:
cc <- binomial(link="cloglog")
d$mu0 <- cc$linkinv(d$eta0)
d$mu <- cc$linkinv(d$eta)

d$y0 <- rbinom(ntot,prob=d$mu0,size=1)
d$y <- rbinom(ntot,prob=d$mu,size=1)

## A. no random effects (vs glm)
g0 <- glmmadmb(y0~x+offset(offset),data=d,
               family="binomial",link="cloglog")

g0A <- glm(y0~x+offset(offset),data=d,
           family=binomial(link="cloglog"))

mlist <- list(glmmadmb0=g0,glm0=g0A)


t(sapply(mlist,coef))
sapply(mlist,logLik)

## B. random effects (vs glmer)
g1 <- glmmadmb(y~x+(1|f)+offset(offset),data=d,
               family="binomial",link="cloglog")
library(lme4)
g1A <- glmer(y~x+(1|f)+offset(offset),data=d,
             family=binomial(link="cloglog"))
c1A <- fixef(g1A)
r1A <- ranef(g1A)
v1A <- VarCorr(g1A)
L1A <- logLik(g1A)
detach("package:lme4")
summary(g1)

c(logLik(g1),L1A)
p0 <- predict(g1)
pb <- predict(g1,type="response")


### GAMMA/LOG LINK
gshape <- 1.5

cc <- Gamma(link="log")
d$mu0 <- cc$linkinv(d$eta0)
d$mu <- cc$linkinv(d$eta)

d$y0 <- rgamma(ntot,shape=gshape,scale=d$mu0/gshape)
d$y <- rgamma(ntot,shape=gshape,scale=d$mu/gshape)

## glmmadmb vs glm, no random effect
g2 <- glmmadmb(y0~x,data=d,family="gamma",link="log")
g2L <- glm(y0~x,data=d,family=Gamma(link="log"))

coef(g2)
coef(g2L)

g3 <- glmmadmb(y~x+(1|f),data=d,family="gamma",link="log")
library(lme4)
## boom!
g3L <- glmer(y~x+(1|f),data=d,family=Gamma(link="log"))
coef(g3)
fixef(g3L)


## POISSON/identity link
dd <- data.frame(y=rpois(20,lambda=10))
g5 <- glmmadmb(y~1,data=dd,
         start=list(fixed=10),
         family="poisson",link="identity",
         verbose=TRUE)
mean(dd)

if (FALSE) {
### GAUSSIAN/IDENTITY LINK
d$mu0 <- d$eta0
d$mu <- d$eta

d$y0 <- rnorm(ntot,d$mu0,sd=1)
d$y <- rnorm(ntot,d$mu,sd=1)

g4 <- lm(y0~x,data=d)
g4B <- glm(y0~x,data=d)
g4C <- glmmadmb(y0~x,data=d,family="gaussian")
g4C <- glmmadmb(y0~1,data=d,family="gaussian")
}
