## tests of a single random grouping variable
library(glmmADMB)

## random-intercept model
set.seed(101)
nblock <- 10
nrep <- 10
d <- expand.grid(f=factor(LETTERS[1:nblock]),rep=1:nrep)
N <- nrow(d)
d$x <- runif(N)
u <- rnorm(nblock,sd=1)
beta <- c(1,2)
eta <- model.matrix(~x,data=d) %*% beta + u[as.numeric(d$f)]
d$y <- rpois(N,exp(eta))
##g1 <- glmm.admb(y~x,random=~1,group="f",family="poisson",data=d)
g1 <- glmmadmb(y~x+(1|f),family="poisson",data=d)
## identical 
coef(g1)
g1$stdbeta
logLik(g1)
g1$S
g1$sd_S
summary(g1$U[[1]])

## try with data in global workspace, not data=argument
y <- d$y
x <- d$x
f <- d$f
g1GW <- glmmadmb(y~x+(1|f),family="poisson")
rm(y,x,f)

## random intercepts and slopes
set.seed(101)
nblock <- 10
nrep <- 10
d2 <- expand.grid(f=factor(LETTERS[1:nblock]),rep=1:nrep)
N <- nrow(d)
d2$x <- runif(N)
u_f <- rnorm(nblock,sd=1)
u_fx <- rnorm(nblock,sd=0.5)
beta <- c(1,2)
eta <- model.matrix(~x,data=d2) %*% beta + u[as.numeric(d2$f)]+
  u[as.numeric(d2$f)]*d$x
d2$y <- rpois(N,exp(eta))

##old style:
if (FALSE) {
    ## FIXME: fails on r-forge linux x64 tests
    g2 <- glmmadmb(y~x,random=~x|f,family="poisson",data=d2)
    g2A <- glmmadmb(y~x+(x|f),family="poisson",data=d2,
                    admb.opts=admbControl(noinit=FALSE))

    coef(g2)
    logLik(g2)
    g2$U
    summary(fitted(g2))
    g2$S
    g2$sd_S
}

if (FALSE) {
  glmmadmb(y~x+(x|f),family="poisson",data=d2,
           save.dir="tmp",
           admb.opts=admbControl(run=FALSE))
}

if (FALSE) {
    ## FIXME: test fails on windows with current versions (r1696 lme4, r211 glmmADMB)
    library(lme4)
    g2B <- glmer(y~x+(1|f)+(0+x|f),family="poisson",data=d2)
    stopifnot(all.equal(fixef(g2B),coef(g2),tol=3e-5))
    stopifnot(all.equal(unname(unlist(VarCorr(g2B))),unname(diag(g2$S[[1]])),tol=3e-3))
}



g2D <- glmmadmb(y~x+(x|f),family="poisson",data=d2,
               admb.opts=admbControl(maxph=NA))

g2E <- glmmadmb(y~x+(x|f),family="poisson",data=d2,
               admb.opts=admbControl(noinit=FALSE))

if (FALSE) {
  ## this one still doesn't work
  g2C <- glmmadmb(y~x+(x|f),family="poisson",data=d2)
}


