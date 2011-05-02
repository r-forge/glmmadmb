library(glmmADMB)  ## testing version

set.seed(101)
## nblock <- 10
## nrep <- 10

## smaller (5 x 50) for quicker testing
simfun <- function(seed=101,
                   nblock=5,
                   nrep=50,
                   rsd=c(f=1,g=2),
                   beta=c(1,2)) {
  if (!is.null(seed)) set.seed(seed)
  d <- expand.grid(f=factor(LETTERS[1:nblock]),
                   g=factor(letters[1:nblock]),
                   rep=1:nrep)
  N <- nrow(d)
  d$x <- runif(N)
  u_f <- rnorm(nblock,sd=rsd["f"])
  u_g <- rnorm(nblock,sd=rsd["g"])
  eta <- model.matrix(~x,data=d) %*% beta +
    u_f[as.numeric(d$f)]+u_g[as.numeric(d$g)]
  d$y <- rpois(N,exp(eta))
  attr(d,"reff") <- rbind(data.frame(eff="f",block=levels(d$f),u=u_f),
                          data.frame(eff="g",block=levels(d$g),u=u_g))
  d
}
                          
                               
library(lme4)
## kluge, for passing tests until I can get this sorted out
setMethod("VarCorr", signature(x="glmm.admb"), glmmADMB:::VarCorr.glmm.admb)
setMethod("VarCorr", signature(x="summary.glmm.admb"), glmmADMB:::VarCorr.glmm.admb)

d1 <- simfun()
m1 <- glmer(y~x+(1|f)+(1|g),family="poisson",data=d1)
res1 <- cbind(attr(d1,"reff"),glmer_est=c(unlist(ranef(m1))))
g1 <- glmm.admb(y~x+(1|f)+(1|g),family="poisson",data=d1)
res1 <- cbind(res1,glmmADMB_est=c(unlist(ranef(g1))))

plot(res1[,"u"],res1[,"glmmADMB_est"],col=as.numeric(res1[,"eff"]))
abline(a=0,b=1)

summary(g1)
coef(g1)
VarCorr(g1)


dd <- data.frame(expand.grid(x=1:5,eff=c("f","g"),
                       method=c("glmmadmb","glmer") ),
           re = c(unlist(g2$U),unlist(ranef(m2))))
dd2 <- recast(dd,...~method,id.var=1:3)
dd2$re <- c(matrix(c(u_f,u_g),nrow=2,byrow=TRUE))
library(ggplot2)
ggplot(dd2,aes(x=re,y=glmer))+geom_point()+facet_grid(.~eff)+
  geom_abline(intercept=0,slope=1)+
  geom_point(aes(y=glmmadmb),colour="red")
##

## random effects are PROPORTIONAL but not IDENTICAL:
##  different scaling in each case?
##  neither matches true values very well
plot(u_f,g2$U$f,xlab="true",ylab="estimated")
abline(a=0,b=1)
plot(g2$U$g,u_g)
abline(a=0,b=1)

###############

d2 <- simfun(nrep=200)
m3 <- glmer(y~x+(1|f)+(1|g),family="poisson",data=d2)
cbind(attr(d2,"reff"),est=c(unlist(ranef(m3))))

g2 <- glmm.admb(y~x+(1|f)+(1|g),family="poisson",data=d)


d3 <- simfun(nrep=500)
m4 <- glmer(y~x+(1|f)+(1|g),family="poisson",data=d3)
ranef(m4)

if (require(R2admb)) {
  ## it was easier for me to do this via R2admb ...
  Zf <- model.matrix(~f-1,data=d1)
  Zg <- model.matrix(~g-1,data=d1)
  X <- model.matrix(~y,data=d1)
  admb_dat <- list(Zf=Zf,Zg=Zg,X=X,y=d1$y,nobs=nrow(d1))
                 
  setup_admb()
  if (!file.exists("g3.RData")) {
  g3 <- do_admb("crossed",data=admb_dat,
                params=list(beta=c(1,0),
                  sigma_f=1,sigma_g=1),
                re=TRUE,
                re_vectors=c(u_f=length(levels(d1$f)),
                  u_g=length(levels(d1$g))),
                checkparam="write",
                checkdata="write",
                ## extra.args="-l1 10000000 -l2 100000000 -l3 10000000 -nl1 10000000",
                extra.args="-ndb 2 -l1 1000000",
                verbose=TRUE)
  save("g3",file="g3.RData")
} else load("g3.RData")

  
}



