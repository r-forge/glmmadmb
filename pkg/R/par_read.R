par_read <- function(name)
{
  filename <- if(tools::file_ext(name)=="par") name else paste(name,".par",sep="")

  tmp <- scan(filename, what="", quiet=TRUE)
  tmp2 <- split(tmp, cumsum(tmp=="#"))
  x <- tmp2[-1]

  for(i in 1:length(x))
  {
    y <- x[[i]]
    n <- nchar(y[2])
    x[[i]] <- as.numeric(y[-(1:2)])
    names(x)[i] <- substring(y[2], 1, n-1)
  }

  x$loglik <- -as.numeric(tmp2[[1]][11])
  x$gradient <- -as.numeric(tmp2[[1]][16])

  return(x)
}
