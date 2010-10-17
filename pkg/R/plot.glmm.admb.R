plot.glmm.admb <- function(x, ...)
{
  plot(x$fitted, x$residuals, xlab="Fitted values", ylab="Residuals", ...)
}
