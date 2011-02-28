library(glmmADMB)  ## testing version

set.seed(101)
## nblock <- 10
## nrep <- 10

## smaller for quicker testing
nblock <- 5
nrep <- 50
d <- expand.grid(f=factor(LETTERS[1:nblock]),
                 g=factor(letters[1:nblock]),
                 rep=1:nrep)
N <- nrow(d)
d$x <- runif(N)
u_f <- rnorm(nblock,sd=1)
u_g <- rnorm(nblock,sd=2)
beta <- c(1,2)
eta <- model.matrix(~x,data=d) %*% beta + u_f[as.numeric(d$f)]+u_g[as.numeric(d$g)]
d$y <- rpois(N,exp(eta))

g2 <- glmm.admb(y~x+(1|f)+(1|g),family="poisson",data=d)

summary(g2)
coef(g2)
VarCorr(g2)

## 
plot(u_f,g2$U$f,xlab="true",ylab="estimated")
abline(a=0,b=1)
plot(g2$U$g,u_g)
abline(a=0,b=1)
