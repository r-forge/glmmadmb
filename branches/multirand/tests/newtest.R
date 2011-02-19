(fm <- glmm.admb(y~Base*trt+Age+Visit+(Visit|subject),
                 data=epil2, family="nbinom",verbose=TRUE))


## Dobson (1990) Page 93: Randomized Controlled Trial :
d <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                outcome=gl(3,1,9),
                treatment <- gl(3,3))
glm.D93 <- glm(counts ~ outcome + treatment,
               data=d,
               family=poisson())

glmm.admb(counts~outcome+treatment,family="poisson",data=d)


## single random effect
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

g1 <- glmm.admb(y~x+(1|f),family="poisson",data=d)


set.seed(101)
nblock <- 10
nrep <- 10
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

