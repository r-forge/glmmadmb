library(glmmADMB)
if (is.null(testLevel <- Sys.getenv("GLMMADMB_TEST_LEVEL"))) testLevel <- 1
## slow, only run when requested
if (testLevel>1) {
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  logBroodSize=log(BroodSize),
                  NCalls=SiblingNegotiation)
fit_zipoiss <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                        offset(logBroodSize)+(1|Nest),
                        data=Owls,
                        zeroInflation=TRUE,
                        family="poisson")
sessionInfo()
}
