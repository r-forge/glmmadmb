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
coef(g1)
g1$stdbeta
logLik(g1)
g1$S
g1$sd_S
summary(g1$U$f) ## not quite identical

library(lme4)
m1 <- glmer(y~x+(1|f),family="poisson",data=d)
r1 <- ranef(g1)[[1]]
r2 <- ranef(m1)[[1]][[1]]
## ranef(mer) returns DATA FRAME so need additional unwrapping
plot(r1,r2,xlab="glmmADMB random effects",ylab="glmer random effects")
abline(a=0,b=1)
## proportional but 'shrunk' slightly differently
summary(lm(r2~r1)) ## slope=0.9387

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

g2 <- glmmadmb(y~x+(x|f),family="poisson",data=d2)
coef(g2)
logLik(g2)
g2$U$f ## not quite identical
summary(fitted(g2))
g2$S
g2$sd_S
summary(g2)

epil2$subject <- factor(epil2$subject)
fm <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                 data=epil2, family="nbinom")
