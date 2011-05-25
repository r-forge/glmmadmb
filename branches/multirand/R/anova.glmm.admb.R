anova.glmmadmb <- function(object, ...)
{
  objects <- list(object, ...)

  if(length(objects) < 2)
    stop("Two or more model fits required.")
  if(length(unique(paste(lapply(objects,function(x) x$random)))) > 1)
    stop("Random effects are not identical")

  npar <- as.numeric(lapply(objects, function(x) x$npar))
  logLik <- as.numeric(lapply(objects, function(x) x$loglik))
  df <- c(NA, diff(npar))
  n2logQ <- 2 * c(NA,diff(logLik))
  P.value <- c(NA, 1-pchisq(n2logQ[-1],df[-1]))
  table <- data.frame(npar, logLik, df, n2logQ, P.value)
  variables <- lapply(objects, function(x) x$fixed)

  dimnames(table) <- list(1:length(objects), c("NoPar","LogLik","Df","-2logQ","P.value"))
  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1:length(objects)), ": ", variables, sep="", collapse="\n")

  output <- structure(table, heading=c(title,topnote), class=c("anova","data.frame"))

  return(output)
}
