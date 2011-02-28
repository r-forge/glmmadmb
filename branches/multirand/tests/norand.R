library(glmmADMB)

## Dobson (1990) Page 93: Randomized Controlled Trial :
d <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                outcome=gl(3,1,9),
                treatment <- gl(3,3))
glm.D93 <- glm(counts ~ outcome + treatment,
               data=d,
               family=poisson())

## FAILS ("error division by zero in solve(dvar_matrix)")
glmm.admb(counts~outcome+treatment,family="poisson",data=d)


