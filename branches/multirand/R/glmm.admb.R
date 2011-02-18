glmm.admb <- function(formula, data, family="poisson", link, corStruct="diag", impSamp=0, easyFlag=TRUE,
                      zeroInflation=FALSE, imaxfn=10, save.dir=NULL, verbose=FALSE)
{
  ## BMB: make this an R temp directory? begin with a .?
  dirname <- if(is.null(save.dir)) "_glmm_ADMB_temp_dir_" else save.dir
  if(!file_test("-d",dirname))
    dir.create(dirname)
  owd <- setwd(dirname); on.exit(setwd(owd))
  if(is.null(save.dir))
    on.exit(unlink(dirname,recursive=TRUE), add=TRUE)

  call <- match.call()
  if(!(corStruct %in% c("diag","full")))
    stop("Argument \"corStruct\" must be either \"diag\" or \"full\"")
  if(!(family %in% c("poisson","nbinom","binomial")))
    stop("Argument \"family\" must be either \"poisson\", \"binomial\"  or \"nbinom\"")

  has_rand <- length(grep("\\|",as.character(formula)[3]))>0
  if (!has_rand && (!missing(impSamp) || !missing(corStruct)))
    stop("No random effects specified: neither \"impSamp\" or \"corStruct\" make sense")
  if(missing(link) && family=="binomial")
    stop("Argument \"link\" must be provided for the \"binomial\" family)")
  ## BMB: make logit link the default?

  ## from glm()
  ## extract x, y, etc from the model formula and frame
  if (missing(data)) 
    data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- get_fixedformula(mf$formula)
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms") # allow model.frame to have updated it

  y <- model.response(mf, "any") # e.g. factors are allowed
  if(length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    if(!is.null(nm)) names(y) <- nm
  }

  n <- if (is.matrix(y)) nrow(y) else length(y)

  offset <- as.vector(model.offset(mf))
  has_offset <- !is.null(offset)
  if (!has_offset) offset <- rep(0,n)

  X <- model.matrix(mt,mf)  ## BMB: allow contrasts?
  p <- ncol(X)
  
  REmat <- process_randformula(formula,data=data)

  m <- sapply(REmat$mmats, function(x) ncol(x))		 # Number of random effects within each crossed term
  M <- length(m)					 # Number of crossed terms
  q <- sapply(REmat$codes, function(x) length(unique(x))) # Number of levels of each crossed term

  # Construct design matrix (Z) for random effects
  for(i in 1:M)
  {
    tmp <- REmat$mmats[[i]]
    colnames(tmp) = paste(names(REmat$mmats)[i],colnames(tmp),sep=":")
    if(i==1)
      Z <- tmp
    else
      Z <- cbind(Z,tmp)
  }

  # Make matrix of pointers into ADMB vector of random effects "u".
  # Rows in II contain pointers for a given data point 
  # First offset the columns of II to point into different segments of u
  II = matrix(rep(q,m),nrow=n,ncol=sum(m),byrow=T)
  if(sum(m)>1)
    for(i in 2:sum(m))
      II[,i] = II[,i] + II[,i-1] 
  II = cbind(0,II[,-ncol(II)])
  ii = 1
  for(i in 1:M)		# Fill in the actual random effects codes
    for(j in 1:m[i])
    {
      II[,ii] = II[,ii] + REmat$codes[[i]]
      ii = ii+1
    }

  # Splits u into "correlation blocks" (no way yet of specifying uncorrelated random effects)
  cor_block_start = cumsum(m)-m[1]+1
  cor_block_stop = cumsum(m)
  numb_cor_params = max(sum(m*(m-1)/2),1)	# For technical reasons we need at least 1 parameter: see glmmadmb.tpl
  
  cmdoptions = paste("-maxfn 500", if(impSamp==0) "" else paste("-is",impSamp))
  dat_list = list(n=n, y=y, p=p, X=X, M=M, q=q, m=m, ncolZ=ncol(Z),Z=Z, II=II, 
		   	cor_flag=rep(0,M),
			cor_block_start=cor_block_start,
		   	cor_block_stop=cor_block_stop,
			numb_cor_params=numb_cor_params,
                   	like_type_flag=as.numeric(family=="poisson"||(family=="binomial"&&link=="logit")),
                   	no_rand_flag=as.numeric(!has_rand),
                   	zi_flag=as.numeric(zeroInflation), 
			intermediate_maxfn=10, 
			has_offset=as.numeric(has_offset), 
			offset=offset)
  ## BMB: pz=0.0001 should be clearly specified, possibly made into a control parameter
  pin_list = list(pz=if(zeroInflation) 0.02 else 0.0001, b=numeric(p), tmpL=0.25+numeric(sum(m)),
    tmpL1=0.0001+numeric(numb_cor_params), logalpha=2.0, u=rep(0,sum(m*q)))


  file_name <- "glmmadmb"

  dat_write(file_name, dat_list)
  pin_write(file_name, pin_list)
  std_file <- paste(file_name, ".std", sep="")
  if(file.exists(std_file))
    file.remove(std_file)

  if(.Platform$OS.type == "windows")
  {
    cmd <- paste("\"",system.file("bin","windows",paste(file_name,".exe",sep=""),package="glmmADMB"), "\"", " ", cmdoptions, sep="")
    shell(cmd, invisible=TRUE)
  } else  {
    if (substr(R.version$os,1,6)=="darwin") {
      ## MacOS detection OK?
      ## http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8497.html
      ## do we really need to copy the binaries over to the temp directory, or can we
      ##   run them in situ?
      file.copy(system.file("bin","macos",file_name,package="glmmADMB"),".")
    } else if (R.version$os=="linux-gnu") {
      file.copy(system.file("bin","linux",file_name,package="glmmADMB"),".")
    } else stop("unknown OS detected")
    Sys.chmod(file_name,mode="0755") ## file.copy strips executable permissions????
    cmd2 <- paste("./", file_name, " ", cmdoptions, sep="")
    sys.result <- system(cmd2,intern=!verbose)
    unlink(file_name)
  }

  if(!file.exists(std_file))
    stop("The function maximizer failed")
  tmp <- read.table(paste(file_name,"std",sep="."), skip=1,as.is=TRUE)
  ## BMB: if 'std dev' were written without a space we could use header=TRUE
  tmpindex <- tmp[,2]
  out <- list(n=n, q=q, fixed=fixed, call=call, group=group, family=family, corStruct=corStruct, impSamp=impSamp,
              easyFlag=easyFlag, zeroInflation=zeroInflation)
  if(zeroInflation)
    out$pz <- as.numeric(tmp[tmpindex=="pz", 3])
  out$b <- as.numeric(tmp[tmpindex=="real_beta", 3])
  out$stdbeta <- as.numeric(tmp[tmpindex=="real_beta", 4])
  names(out$stdbeta) <- names(out$b) <- colnames(X)

  if(family == "nbinom")
  {
    out$alpha <- as.numeric(tmp[tmpindex=="alpha", 3])
    out$sd_alpha <- as.numeric(tmp[tmpindex=="alpha", 3])
  }
  if(has_rand) {
    ## BMB: fixme!
    out$S <- matrix(as.numeric(tmp[tmpindex=="S",3]), nrow=m, dimnames=list(colnames(Z),colnames(Z)))
    out$sd_S <- matrix(as.numeric(tmp[tmpindex=="S",4]), nrow=m, dimnames=list(colnames(Z),colnames(Z)))
    if(corStruct == "diag")
    {
      out$S[below(m,TRUE)|t(below(m,TRUE))] <- 0
      out$sd_S[below(m,TRUE)|t(below(m,TRUE))] <- 0
    }
    ## out$random <- random
  }
  if(!missing(link))
    out$link <- link
  if(has_rand)
  {
    U <- matrix(as.numeric(tmp[tmpindex=="u",3]), ncol=m, byrow=TRUE,
                dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    out$U <- U
    out$sd_U <- matrix(as.numeric(tmp[tmpindex=="u",4]), ncol=m, byrow=TRUE,
                       dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
  }
  else
  {
    U <- matrix(rep(0,q), ncol=1, byrow=TRUE)
  }

  mu <- as.numeric(X %*% out$b)
  ## BMB: doesn't include influence of random effects?
  lambda <- 0
  for(i in 1:n)
    lambda[i] <- exp(mu[i] + sum(Z[i,]*U[II[i],]))
  if(family == "binomial")
    out$fitted <- lambda / (1+lambda)
  else
    out$fitted <- lambda
    out$sd.est <- switch(family,
                  poisson=sqrt(lambda),
                  nbinom=sqrt(lambda*(1+lambda/out$alpha)),
                  binomial=sqrt(out$fitted*(1-out$fitted)))

  out$residuals <- as.numeric(y-lambda)
  tmp <- par_read(file_name)
  out$npar <- tmp$npar ## as.numeric(scan(paste(file_name,".par",sep=""), what="", quiet=TRUE)[6])
  ## BMB: should this be total number of parameters or number of fixed parameters?
  out$loglik <- tmp$loglik
  out$gradloglik <- tmp$gradient
  nfixpar <- length(out$b)
  ## drop cors that don't correspond to fixed-effect parameters
  out$corMat <- tmp$cor[1:nfixpar,1:nfixpar]

  out$conv <- 0
  out$convmsg <- ""
  
  if (abs(out$gradloglik) >= 0.001) {
    out$convmsg <- paste("log-likelihood of gradient=",out$gradloglik)
    out$conv <- 1
    warning("Convergence failed:",out$convmsg)
  }

  ## unnecessary??
  ## out$fitted = out$fitted[backwards_key]
  ## out$residuals = out$residuals[backwards_key]

  class(out) <- "glmm.admb"

  return(out)
}
