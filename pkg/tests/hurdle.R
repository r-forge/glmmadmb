set.seed(1001)

data("bioChemists",package="pscl")
## library(pscl)
## h1 <- hurdle(art~fem+mar+kid5+phd+ment,data=bioChemists)

psclcoef <- structure(c(0.671139312871508,
                        -0.228582657955377, 0.0964849892769541, 
                        -0.142187559329214, -0.0127263725842253,
                        0.0187454797396211, 
                        0.236796012435245, -0.251151128617988,
                        0.326233583613064, -0.285248715785051, 
                        0.0222193970990452, 0.0801213546865395),
                      .Names = c("count_(Intercept)", 
                        "count_femWomen", "count_marMarried",
                        "count_kid5", "count_phd", 
                        "count_ment", "zero_(Intercept)", "zero_femWomen",
                        "zero_marMarried", "zero_kid5",
                        "zero_phd", "zero_ment"))

bb <- subset(bioChemists,art>0)
library(glmmADMB)
g1 <- glmmadmb(art~fem+mar+kid5+phd+ment,
               family="truncpoiss",link="log",data=bb)
cp1 <- psclcoef[grep("count_",names(psclcoef))]
cg1 <- coef(g1)
stopifnot(abs(cp1-cg1)<3e-6)

bioChemists <- transform(bioChemists,nz=as.numeric(art>0))
g2 <- glmmadmb(nz~fem+mar+kid5+phd+ment,
               family="binomial",data=bioChemists)
cp2 <-  psclcoef[grep("zero_",names(psclcoef))]
cg2 <- coef(g2)
stopifnot(abs(cp2-cg2)<2e-5)
