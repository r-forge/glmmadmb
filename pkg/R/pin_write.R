"pin_write" <- function(name, L)
{
  n = nchar(name)
  if(substring(name,n-3,n) == ".pin")
    file_name = name
  else
    file_name = paste(name, ".pin", sep="")
  cat("# \"", name, ".pin\" produced by pin_write() from ADMButils; ", date(), "\n", file=file_name, sep="")
  for(i in 1:length(L))
  {
    x = L[[i]]
    if(data.class(x) == "numeric")
      cat("#", names(L)[i], "\n", L[[i]], "\n\n", file=file_name, append=TRUE)
    if(data.class(x) == "matrix")
    {
      cat("#", names(L)[i], "\n", file=file_name, append=TRUE)
      write.table(L[[i]], , col=FALSE, row=FALSE, quote=FALSE, file=file_name, append=TRUE)
      cat("\n", file=file_name, append=TRUE)
    }
  }
}
