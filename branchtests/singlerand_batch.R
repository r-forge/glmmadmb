## should be run with --vanilla

## simulated data, intercept-only RE
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

## simulated data, slope+intercept RE
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

##
data(epil2,package="glmmADMB")

## model 1:

## old glmmADMB
library(glmmADMB.old)  ## version 0.5-2
t0_old <- system.time(g0_old <- glmm.admb(y~1,random=~1,
                                          group="f",family="poisson",data=d))
t1_old <- system.time(g1_old <- glmm.admb(y~x,random=~1,
                                          group="f",family="poisson",data=d))
t2_old <- system.time(g2_old <- glmm.admb(y~x,random=~x,
                                          group="f",family="poisson",data=d))
t3_old <- system.time(g3_old <- glmm.admb(y~1,random=~1,group="f",
                                          family="poisson",data=d2))
t4_old <- system.time(g4_old <- glmm.admb(y~x,random=~1,group="f",
                                          family="poisson",data=d2))
t5_old <- system.time(g5_old <- glmm.admb(y~x,random=~x,group="f",
                                          family="poisson",data=d2))
t6_old <- system.time(g6_old <- glmm.admb(y~Base*trt+Age+Visit, 
                                          random=~Visit, group="subject",
                                          data=epil2, family="nbinom"))
## epil2 data fitted with Poisson, for comparison to lme4
t7_old <- system.time(g7_old <- glmm.admb(y~Base*trt+Age+Visit, 
                                          random=~Visit, group="subject",
                                          data=epil2, family="poisson"))
save.image("singlerand_batch.RData")

detach("package:glmmADMB.old")

library("glmmADMB")
t0_new <- system.time(g0_new <- glmm.admb(y~1+(1|f),
                                          family="poisson",data=d))
t1_new <- system.time(g1_new <- glmm.admb(y~x+(1|f),
                                          family="poisson",data=d))
t2_new <- system.time(g2_new <- glmm.admb(y~x+(x|f),
                                          family="poisson",data=d))
t3_new <- system.time(g3_new <- glmm.admb(y~1+(1|f),
                                          family="poisson",data=d2))
t4_new <- system.time(g4_new <- glmm.admb(y~x+(1|f),
                                          family="poisson",data=d2))
t5_new <- system.time(g5_new <- glmm.admb(y~x+(x|f),
                                          family="poisson",data=d2))
t6_new <- system.time(g6_new <- try(glmm.admb(y~Base*trt+Age+Visit+
                                          (Visit|subject),
                                          data=epil2, family="nbinom")))
## epil2 data fitted with Poisson, for comparison to lme4
t7_new <- system.time(g7_new <- glmm.admb(y~Base*trt+Age+Visit+ 
                                          (Visit|subject),
                                          data=epil2, family="poisson"))
save.image("singlerand_batch.RData")
detach("package:glmmADMB")

library("lme4")
t0_lme4 <- system.time(g0_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d))
t1_lme4 <- system.time(g1_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d))
t2_lme4 <- system.time(g2_lme4 <- glmer(y~x+(x|f),
                                          family="poisson",data=d))
t3_lme4 <- system.time(g3_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d2))
t4_lme4 <- system.time(g4_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d2))
t5_lme4 <- system.time(g5_lme4 <- glmer(y~x+(x|f),
                                          family="poisson",data=d2))
## t6_lme4 <- system.time(g6_lme4 <- try(glmer(y~Base*trt+Age+Visit+
##                                           (Visit|subject),
##                                           data=epil2, family="nbinom")))
## epil2 data fitted with Poisson, for comparison to lme4
t7_lme4 <- system.time(g7_lme4 <- glmer(y~Base*trt+Age+Visit+ 
                                          (Visit|subject),
                                          data=epil2, family="poisson"))
save.image("singlerand_batch.RData")


