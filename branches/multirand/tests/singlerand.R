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

g1 <- glmmadmb(y~x+(1|f),family="poisson",data=d)
g1B <- glmmadmb(y~x,random=~(1|f),family="poisson",data=d)
g1C <- glmmadmb(y~x,random=~1|f,family="poisson",data=d)
coef(g1)
fixef(g1)
summary(g1)
logLik(g1)
ranef(g1)
AIC(g1)

library(lme4)
m1 <- glmer(y~x+(1|f),family="poisson",data=d)
r1 <- unname(unlist(ranef(g1)))
r2 <- unname(unlist(ranef(m1)))
## ranef() returns DATA FRAME so need additional unwrapping
plot(r1,r2,xlab="glmmADMB random effects",ylab="glmer random effects")
abline(a=0,b=1)

stopifnot(all.equal(fixef(g1),fixef(m1),tolerance=1e-6))
stopifnot(all.equal(r1,r2,tolerance=2e-4))

## random intercepts and slopes
set.seed(101)
nblock <- 10
nrep <- 10
d2 <- expand.grid(f=factor(LETTERS[1:nblock]),rep=1:nrep)
N <- nrow(d2)
d2$x <- runif(N)
u_f <- rnorm(nblock,sd=1)
u_fx <- rnorm(nblock,sd=0.5)
beta <- c(1,2)
eta <- model.matrix(~x,data=d2) %*% beta + u[as.numeric(d2$f)]+
  u[as.numeric(d2$f)]*d2$x
d2$y <- rpois(N,exp(eta))

(g2 <- glmmadmb(y~x+(x|f),family="poisson",data=d2))
summary(g2)
m2 <- glmer(y~x+(1|f)+(0+x|f),family="poisson",data=d2)
r1 <- unname(unlist(ranef(g2)))
r2 <- unname(unlist(ranef(m2)))

stopifnot(all.equal(fixef(g2),fixef(m2),tolerance=1e-4))
stopifnot(all.equal(r1,r2,tolerance=1.5e-3))

epil2$subject <- factor(epil2$subject)
(fm <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                 data=epil2, family="nbinom"))
summary(fm)
