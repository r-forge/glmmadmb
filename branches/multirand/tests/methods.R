library(glmmADMB)
data(OwlModel)
summary(OwlModel)
coef(OwlModel)
stdEr(OwlModel)
vcov(OwlModel)
AIC(OwlModel)
## AICc: implemented in MuMIn, AICcmodavg, bbmle
## glmmADMB objects work with MuMIn and bbmle versions, not
##   with AICcmodavg
## FIXME: write to AICcmodavg maintainer?
