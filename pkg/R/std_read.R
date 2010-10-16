"std_read" <- function(name)
{
  n = nchar(name)

  if(substring(name,n-3,n) == ".std")
    file_name = name
  else
    file_name = paste(name, ".std", sep="")

  tmp = read.table(file_name, skip=1)
  est = tmp[,3]
  names(est) = tmp[,2]
  std = tmp[,3]
  names(std) = tmp[,2]
  tmp = scan(paste(name,".par",sep=""), what="", quiet=TRUE)
  loglik = as.numeric(tmp[11])
  grad = as.numeric(tmp[16])

  list(est, std, loglik, grad)
}
