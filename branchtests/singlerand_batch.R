## should be run with --vanilla

## compare all known examples with a SINGLE random variable that
## glmm.admb (and other models) are capable of.

## simulated data, intercept-only RE
simfun <- function(seed=101,nblock=10,nrep=10,sd=c(int=1,slope=0),beta=c(1,2)) {
  if (!is.null(seed)) set.seed(seed)
  d <- expand.grid(f=factor(LETTERS[1:nblock]),rep=1:nrep)
  N <- nrow(d)
  d$x <- runif(N)
  u_f1 <- rnorm(nblock,sd=sd["int"])
  u_fx <- rnorm(nblock,sd=sd["slope"])
  eta <- model.matrix(~x,data=d) %*% beta + u_f1[as.numeric(d$f)]+
    u_fx[as.numeric(d$f)]*d$x    
  d$y <- rpois(N,exp(eta))
  d
}

## simulated data, slope+intercept RE
d <- simfun()
d2 <- simfun(sd=c(int=1,slope=0.5))

cdata <- read.csv("culcitalogreg.csv",
              colClasses=c(rep("factor",2),
                "numeric",
                rep("factor",6)))
cdata$block <- factor(cdata$block,levels=1:10)

cdata2 <- cdata[order(cdata$block),]

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

t8_old <- system.time(g8_old <- glmm.admb(predation~ttt,
                                          random=~1,
                                          group="block",family="binomial",
                                          link="logit",
                                          data=cdata,
                                          save.dir="g8oldtest"))

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
if (FALSE) {
  ## SKIP: binomial not yet implemented in new glmmADMB!
t8_new <- system.time(g8_new <- glmm.admb(predation~ttt+(1|block),
                                          family="binomial",data=cdata,
                                          save.dir="g8newtest"))
t8_new2 <- system.time(g8_new2 <- glmm.admb(predation~ttt+(1|block),
                                            family="binomial",data=cdata2,
                                            save.dir="g8newtest2"))
}

library("lme4")
t0_lme4 <- system.time(g0_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d))
t1_lme4 <- system.time(g1_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d))
t2_lme4 <- system.time(g2_lme4 <- glmer(y~x+(1|f)+(0+x|f),
                                          family="poisson",data=d))
t3_lme4 <- system.time(g3_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d2))
t4_lme4 <- system.time(g4_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d2))
t5_lme4 <- system.time(g5_lme4 <- glmer(y~x+(1|f)+(0+x|f),
                                          family="poisson",data=d2))
## t6_lme4 <- system.time(g6_lme4 <- try(glmer(y~Base*trt+Age+Visit+
##                                           (Visit|subject),
##                                           data=epil2, family="nbinom")))
## epil2 data fitted with Poisson, for comparison to lme4
t7_lme4 <- system.time(g7_lme4 <- glmer(y~Base*trt+Age+Visit+ 
                                          (1|subject)+(0+Visit|subject),
                                          data=epil2, family="poisson"))

t8_lme4 <- system.time(g8_lme4 <- glmer(predation~ttt+(1|block),family=binomial,data=cdata))

save.image("singlerand_batch.RData")
detach("package:glmmADMB")

library("lme4")
t0_lme4 <- system.time(g0_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d))
t1_lme4 <- system.time(g1_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d))
t2_lme4 <- system.time(g2_lme4 <- glmer(y~x+(1|f)+(0+x|f),
                                          family="poisson",data=d))
t3_lme4 <- system.time(g3_lme4 <- glmer(y~1+(1|f),
                                          family="poisson",data=d2))
t4_lme4 <- system.time(g4_lme4 <- glmer(y~x+(1|f),
                                          family="poisson",data=d2))
t5_lme4 <- system.time(g5_lme4 <- glmer(y~x+(1|f)+(0+x|f),
                                          family="poisson",data=d2))
## t6_lme4 <- system.time(g6_lme4 <- try(glmer(y~Base*trt+Age+Visit+
##                                           (Visit|subject),
##                                           data=epil2, family="nbinom")))
## epil2 data fitted with Poisson, for comparison to lme4
t7_lme4 <- system.time(g7_lme4 <- glmer(y~Base*trt+Age+Visit+ 
                                          (1|subject)+(0+Visit|subject),
                                          data=epil2, family="poisson"))

t8_lme4 <- system.time(g8_lme4 <- glmer(predation~ttt+(1|block),family=binomial,data=cdata))

save.image("singlerand_batch.RData")

## only does Poisson/binomial random-intercept models (2,5,6,7 not possible)
library(glmmML)
t0_glmmML <- system.time(g0_glmmML <- glmmML(y~1,cluster=f,
                                             family="poisson",data=d,
                                             control=list(maxit=1000)))
t1_glmmML <- system.time(g1_glmmML <- glmmML(y~x,cluster=f,
                                             family="poisson",data=d,
                                             control=list(maxit=1000)))
t3_glmmML <- system.time(g3_glmmML <- glmmML(y~1,cluster=f,
                                          family="poisson",data=d2,
                                             control=list(maxit=1000)))
t4_glmmML <- system.time(g4_glmmML <- glmmML(y~x,cluster=f,
                                          family="poisson",data=d2))

t8_glmmML <- system.time(g8_glmmML  <- glmmML(predation~ttt,
                                              family=binomial,
                                              data=cdata,cluster=block))
## convergence problem on models 1, 3
detach("package:glmmML")

library(lme4a)
t0_lme4a <- system.time(g0_lme4a <- glmer(y~1+(1|f),
                                          family="poisson",data=d))

save.image("singlerand_batch.RData")
