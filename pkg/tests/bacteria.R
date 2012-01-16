data(bacteria,package="MASS")

bacteria$present <- as.numeric(bacteria$y)-1

library(lme4)
gfit <- glmer(present ~ trt + I(week > 2)+(1|ID),family = "binomial", data = bacteria)

library(glmmADMB)
bfit <-  glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria)
stopifnot(all.equal(bfit$phi,
                    glmmADMB:::make_phi(model.matrix(~trt+I(week>2),data=bacteria)),
                                        tol=1e-4))


gsd <- attr(lme4::VarCorr(gfit)[[1]],"stddev")
f <- lme4::fixef(gfit)
u <- unlist(lme4::ranef(gfit)[[1]])/gsd
bfit2 <-glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
                 start =list(fixed=f,RE_sd=log(gsd),u=u),
                 verbose=TRUE)
slist <- list(fixed=f,RE_sd=log(gsd),u=u)
bfit3 <-glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
                 start=slist,
                 extra.args="-phase 5",
                 verbose=TRUE)


## 9.6130683
