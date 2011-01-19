library(glmmADMB)
library(lme4)

## started out as an attempt to test residuals calculations
## between lme4 and

## would like to do a simpler example (e.g. glm only, from ?glm)
## but glmm.admb requires a random effect

## zero-inflated negative binomial is probably best for fitting data:
## glmer can do Poisson-lognormal, not NB, and not zero-inflation

## ... so fit Poisson with RE, even though it's not a great fit to the data

data(Owls)
Owls$Lbroodsize <- log(Owls$BroodSize)

OwlModel_poiss.glmer <- glmer(SiblingNegotiation ~ FoodTreatment * SexParent +
                              (1|Nest)+offset(Lbroodsize),
                              data=Owls, family=poisson)


OwlModel_poiss.admb <- glmm.admb(SiblingNegotiation ~ FoodTreatment * SexParent,
                                 random=~1,
                                 group="Nest",
                                 offset="Lbroodsize",
                                 data=Owls, family="poisson",
                                 easyFlag=FALSE,
                                 verbose=TRUE)

OwlModel_poiss.glmer2 <- glmer(SiblingNegotiation ~ FoodTreatment * SexParent +
                              (1|Nest)+offset(Lbroodsize),
                              data=Owls, family=poisson,
                               start=list(fixed=coef(OwlModel_poiss.admb)))

## "Estimated covariance matrix may not be positive definite"
logLik(OwlModel_poiss.admb)
logLik(OwlModel_poiss.glmer)

if (require(ggplot2)) {
  ca <- fixef(OwlModel_poiss.glmer)
  cg <- coef(OwlModel_poiss.admb)
  d <- rbind(data.frame(par=names(ca),est=ca,pkg="lme4"),
             data.frame(par=names(cg),est=cg,pkg="glmmADMB"))
  levels(d$par) <- c("food","food:sex","intercept","sex")
  qplot(est,par,data=d,colour=pkg)
}

r_glmer <- residuals(OwlModel_poiss.glmer)
r_admb <- residuals(OwlModel_poiss.admb)

plot(r_glmer,r_admb)                       

