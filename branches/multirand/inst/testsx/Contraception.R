library(mlmRev)
library(glmmADMB)
library(lme4)

fm1 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
fm2 <- glmer(use ~ urban+age+livch+(urban|district), Contraception, binomial)

gm1 <- glmm.admb(use ~ urban+age+livch+(1|district), Contraception, "binomial",verbose=TRUE)
gm2 <- glmm.admb(use ~ urban+age+livch+(urban|district), Contraception, "binomial")
