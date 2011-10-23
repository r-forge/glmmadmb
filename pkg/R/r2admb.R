## everything below here copied from R2admb
## modified from R2admb: propagate back?
read_psv <- function(f,names) {
  f <- tolower(f) ## arghv
  fn <- paste(f,"psv",sep=".")
  if (!file.exists(fn)) stop("no PSV file found")
  ans <- read_admbbin(fn)
  if (missing(names)) names <- paste("V",seq(ncol(ans)),sep="")
  colnames(ans) <- names
  ans <- as.data.frame(ans)
  ans
}

read_chunk <- function(fn,sep="^#",maxlines=1000) {
  end <- FALSE
  ans <- character(maxlines)
  i <- 1
  has_sep <- function(x) length(grep(sep,x))>0
  while (!end) {
    tmp <- readLines(fn,n=1)
    if (i>1 && has_sep(tmp)) {
      end=TRUE
      pushBack(tmp,fn)
    } else if (length(tmp)==0) {
      end=TRUE
    } else {
      ans[i] <- tmp
      i <- i+1
    }
  }
  ans[1:(i-1)]
}

read_hst <- function(fn) {
  fn <- paste(fn,"hst",sep=".")
  if (!file.exists(fn)) {
    warning("file ",fn," not found: returning NULL")
    return(NULL)
  }
  f <- file(fn,open="r")
  r <- list()
  repeat {
    chunk <- read_chunk(f)
    ## cat(length(chunk),":",chunk[1],"\n") 
    if (length(chunk)<2) break
    r <- c(r,list(chunk))
  }
  labs <- sapply(r,"[[",1)
  ## single values
  w <- c(1:2,8,10)
  r[w] <- lapply(r[w],
                 function(x) as.numeric(x[2]))
  names(r)[w] <- c("sampsize","stepsize_scale","npars","rseed")
  ## vectors of value
  w <- c(3:7,9)
  r[w] <- lapply(r[w],function(x) 
                 as.numeric(strsplit(gsub("^ +","",x[2])," ")[[1]]))
  names(r)[w] <- c("stepsizes","means","sdevs","lower","upper","mcmcparms")
  ## r$npars is NOT RELIABLE! use length(stepsizes instead)
  ## parameter matrices
  w <- 11:(10+length(r$stepsizes))
  r[w] <- lapply(r[w],
                 function(z) {
                   do.call(rbind,
                           lapply(z[-c(1,length(z))],
                                  function(x) {as.numeric(strsplit(x," ")[[1]])}))})
  names(r)[w] <- gsub("^#","",
                      gsub("\\[([0-9]+)\\]",".\\1",
                           gsub("; *//.*","",labs[w])))
  ans <- c(r[1:10],hists=list(r[w]))
  class(ans) <- "admb_hist"
  ans
}

## read a "standard" ADMB format binary file into R:
##  standard format is: 1 integer describing number
##  of (double) values per vector
read_admbbin <- function(fn) {
  f <- file(fn,open="rb")
  nv <- readBin(f,"int")
  fs <- file.info(fn)$size
  isize <- 4; dsize <- 8
  m <- matrix(readBin(f,"double",n=(fs-isize)/dsize),byrow=TRUE,
         ncol=nv)
  m
}
         
