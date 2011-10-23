d <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                outcome=gl(3,1,9),
                treatment <- gl(3,3))

library(glmmADMB)
g2 <- glmmadmb(counts~outcome+treatment,family="poisson",data=d)
g2 <- glmmadmb(counts~outcome+treatment,family="poisson",data=d,
               save.dir="dirtst")
v <- unlink("dirtst",recursive=TRUE)
stopifnot(v==0)
