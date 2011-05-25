set.seed(1001)
nb <- c(5,8,10)
d <- data.frame(x=runif(prod(nb)),expand.grid(LETTERS[1:nb[1]],letters[1:nb[2]]))
d$inter <- interaction(d$Var1,d$Var2)
sdvec <- c(2,1)
u1 <- rnorm(nb[1],sd=sdvec[1])
u2 <- rnorm(nb[1]*nb[2],sd=sdvec[2])
d$eta <- with(d,1+0.5*x+u1[Var1]+u2[inter])
d$y <- rpois(prod(nb),exp(d$eta))

names(d)[c(2,4)] <-  c("f1","f2")
d <- subset(d,select=c(x,f1,f2,y))

library(lme4)
t1 <- system.time(g1 <- glmer(y~x+(1|f1/f2),data=d,family=poisson))

library(glmmADMB)
t2 <- system.time(g2 <- glmm.admb(y~x+(1|f1/f2),data=d,family="poisson"))
