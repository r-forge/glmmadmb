
## code for extracting random effects
## note: for this to work we shouldn't have anything funky
##  with nested parentheses.  Also don't know what will
##  happen if fixed effects are strewn between random effects:

## multiple terms of the form
##  (<expr1>|<expr2>)
## where <expr1> describes the covariates affected by/interacting with
## the random effect and <expr2> describes the structure of the random effect

get_fixedformula <- function(f) {
  lchar <- as.character(f[2])
  rchar <- as.character(f[3]) ## RHS
  offsetstr <- ""
  if (length(grep("offset\\(",rchar))>0) {
    ## protect/remove offset
    offsetstr <- gsub(".*(\\+ *offset\\([^)]+\\)).*","\\1",rchar)
    rchar <- gsub("offset\\([^)]+\\)","",rchar)
  } 
  rchar <- gsub("\\([^)]+\\)","",rchar) ## parentheses
  rchar <- gsub("(\\+ *\\+ *)+","+",rchar) ## duplicated +
  rchar <- gsub(" *\\+ *$","",rchar) ## terminating +
  as.formula(paste(lchar,"~",rchar,offsetstr))
}

process_randformula <- function(f,data) {
  rchar <- as.character(f[3]) ## RHS
  ## drop bits before first/after last parenthesis
  rchar <- gsub("^[^()]*\\(","",
                gsub(")[^()]*$","",rchar))
  randbits <- grep("\\|",strsplit(rchar,"[()]")[[1]],value=TRUE)
  splitbits <- strsplit(randbits,"\\|")
                              
  cfun <- function(lbit,mdata) {
    m <- model.matrix(as.formula(paste("~",lbit)),mdata)
    m
  }

  ## here we want to expand RHS and provide a list of indices
  ## into the appropriate factor:
  rfun <- function(rbit,rdata) {
    f <- as.formula(paste("~",rbit,"-1"))
    t <- terms(f,data=rdata)
    ## ugly: "If the answer is parse() you should usually rethink the question" but ??
    labs <- attr(t,"term.labels")
    sapply(labs,
           function(lab) as.numeric(with(data,eval(parse(text=lab)))))
  }
  termnames <- gsub("\\|","_bar_",
                    gsub(":","_int_",
                         gsub("/","_nest_",
                              gsub("\\*","_cross",
                                   gsub(" ","",randbits)))))
  groups <- gsub("^ +","",lapply(splitbits,"[",2))

  nonfactors <- groups[!sapply(data[groups],inherits,"factor")]
  if (length(nonfactors)>0)
    stop("all grouping variables must be factors")
  ## FIXME: say which ones
  
  LHS <- gsub("^ +","",lapply(splitbits,"[",1))
  
  L <- list(mmats=lapply(LHS,cfun,mdata=data),
       codes=lapply(groups,rfun,rdata=data))
  for (i in seq_along(L$mmats)) {
    attr(L$mmats[[i]],"levels") <- levels(with(data,get(groups[i])))
  }
  names(L$mmats) <- groups
  names(L$codes) <- termnames
  L
}

## print tp file (not necessarily 
write_randformula <- function(x,name) {
  require(glmmADMB)
  fn <- if(substring(name,nchar(name)-3)==".dat") {
    name
  } else paste(name,".dat",sep="")
  cat("### design matrices for random effects:\n",file=fn)
  dat_write(name,x$mmat,append=TRUE)
  cat("### indices for random effects:\n",file=fn,append=TRUE)
  dat_write(name,x$codes,append=TRUE)
}

