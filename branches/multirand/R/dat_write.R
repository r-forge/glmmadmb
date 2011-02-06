dat_write <- function(name, L, append=FALSE)
{
  filename <- if(substring(name,nchar(name)-3)==".dat") name else paste(name,".dat",sep="")
  ## filename <- if(tools::file_ext(name)=="dat") name else paste(name,".dat",sep="") # R 2.11 and later

  cat("# \"", name, ".dat\" produced by dat_write() from ADMButils; ",
      date(), "\n", file=filename, sep="", append=append)

  for(i in 1:length(L))
  {
    x <- L[[i]]
    if(data.class(x) == "numeric")
    {
      cat("#", names(L)[i], "\n", L[[i]], "\n\n", file=filename, append=TRUE)
    }
    if(data.class(x) == "matrix")
    {
      cat("#", names(L)[i], "\n", file=filename, append=TRUE)
      write.table(L[[i]], col=FALSE, row=FALSE, quote=FALSE, file=filename, append=TRUE)
      cat("\n", file=filename, append=TRUE)
    }
  }

  invisible(NULL)
}
