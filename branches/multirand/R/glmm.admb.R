mcmc.control <- function(mcmc=1000,
                         mcsave,
                         mcnoscale=FALSE,
                         mcgrope=FALSE,
                         mcmult=1) {
  if (missing(mcsave)) mcsave <- pmax(1,floor(1000/mcmc))
  list(mcmc=mcmc,mcsave=mcsave,mcnoscale=mcnoscale,mcgrope=mcgrope,mcmult=mcmult)
}

mcmc.args <- function(L) {
  argstr <- mapply(function(n,val) {
    if (is.numeric(val)) paste("-",n," ",val,sep="") else
    if (isTRUE(val)) paste("-",val,sep="")
  },names(L),L)
  paste(unlist(argstr),collapse=" ")
}

glmm.admb <- function(formula, data, family="poisson", link,
                      corStruct="diag", impSamp=0, easyFlag=TRUE,
                      zeroInflation=FALSE, ZI_kluge=FALSE,
                      imaxfn=10,
                      mcmc=FALSE,
                      mcmc.opts=mcmc.control(),
                      save.dir=NULL, verbose=FALSE)
{
  ## FIXME: removed commented code after checking
  ## FIXME: make this an R temp directory? begin with a . for invisibility?
  dirname <- if(is.null(save.dir)) "_glmm_ADMB_temp_dir_" else save.dir
  if(!file_test("-d",dirname))
    dir.create(dirname)
  owd <- setwd(dirname); on.exit(setwd(owd))
  if(is.null(save.dir))
    on.exit(unlink(dirname,recursive=TRUE), add=TRUE)

  call <- match.call()
  ## FIXME: corStruct could be a vector, corresponding to the random
  ##    effects in order?
  if(!(corStruct %in% c("diag","full")))
    stop("Argument \"corStruct\" must be either \"diag\" or \"full\"")
  has_rand <- length(grep("\\|",as.character(formula)[3]))>0
  if (!has_rand && (!missing(impSamp) || !missing(corStruct)))
    stop("No random effects specified: neither \"impSamp\" or \"corStruct\" make sense")

  like_type_flag <- switch(family,nbinom=0,poisson=1,binomial=2,gamma=3,
                           stop("unknown family"))

  if (missing(link)) {
    link <- switch(family, binomial="logit", nbinom=, poisson=, gamma="log")
  }
  linkfun <- switch(link,log=log,logit=qlogis,probit=qnorm,
                    stop("unknown link function"))
  ilinkfun <- switch(link,log=exp,logit=plogis,probit=pnorm)

  link_type_flag <- switch(link,log=0,logit=1,probit=2)
  
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

  y <- as.matrix(model.response(mf, "any")) # e.g. factors are allowed
  n <- nrow(y)
  p_y <- ncol(y)

  offset <- as.vector(model.offset(mf))
  has_offset <- !is.null(offset)
  if (!has_offset) offset <- rep(0,n)

  X <- model.matrix(mt,mf)  ## BMB: allow contrasts?
  p <- ncol(X)

  if (has_rand) {
    REmat <- process_randformula(formula,data=data)
    ## Number of random effects within each crossed term
    m <- sapply(REmat$mmats, ncol)
    ## Number of crossed terms
    M <- length(m)
    ## Number of levels of each crossed term
    q <- sapply(REmat$codes, function(x) length(unique(x)))
    names(q) <- names(REmat$mmats)

    ## Construct design matrix (Z) for random effects

    Z <- do.call(cbind,REmat$mmats)
    colnames(Z) <- paste(rep(names(REmat$mmat),m),
                         sapply(REmat$mmat,colnames),sep=":")

    ## Make matrix of pointers into ADMB vector of random effects "u".
    ## Rows in II contain pointers for a given data point 
    ## First offset the columns of II to point into different segments of u
    II <- matrix(rep(q,m),nrow=n,ncol=sum(m),byrow=TRUE)
    if(sum(m)>1)
    for(i in 2:sum(m))
      II[,i] = II[,i] + II[,i-1]
    II <- cbind(0,II[,-ncol(II)])
    ##    tmpf <- function(x) c(0,cumsum(x)[-length(x)])
    ##    II <- t(apply(II,1,tmpf))
    ii <- 1
    for(i in 1:M)		# Fill in the actual random effects codes
      for(j in 1:m[i])
        {
          II[,ii] = II[,ii] + REmat$codes[[i]]
          ii = ii+1
        }
    
    ## Splits u into "correlation blocks" (no way yet of specifying uncorrelated random effects)
    cor_block_start = cumsum(m)-m[1]+1
    cor_block_stop = cumsum(m)
    numb_cor_params = max(sum(m*(m-1)/2),1)	# For technical reasons we need at least 1 parameter: see glmmadmb.tpl
    
  } else {
    ## no random factor: fill in with dummies
    m <- M <- q <- 1
    Z <- matrix(0,ncol=1,nrow=n)
    II <- matrix(rep(q,m),nrow=n,ncol=sum(m),byrow=TRUE)
    cor_block_start <- cor_block_stop <- 1
    numb_cor_params <- 1
  }
  cmdoptions = "-maxfn 500"
  if(impSamp>0) cmdoptions <- paste(cmdoptions,"-is",impSamp)
  if (mcmc) cmdoptions <- paste(cmdoptions,mcmc.args(mcmc.opts))
  dat_list = list(n=n, p_y=p_y,
    y=y, p=p, X=X, M=M, q=q, m=m, ncolZ=ncol(Z),Z=Z, II=II, 
    cor_flag=rep(0,M),
    cor_block_start=cor_block_start,
    cor_block_stop=cor_block_stop,
    numb_cor_params=numb_cor_params,
    like_type_flag=like_type_flag,
    link_type_flag=link_type_flag,
    ## as.numeric(family=="poisson"||(family=="binomial"&&link=="logit")),
    no_rand_flag=as.numeric(!has_rand),
    zi_flag=as.numeric(zeroInflation),
    zi_kluge=as.numeric(ZI_kluge),
    intermediate_maxfn=10, 
    has_offset=as.numeric(has_offset), 
    offset=offset)
  ## BMB: pz=0.0001 should be clearly specified, possibly made into a control parameter
  pin_list = list(pz=if(zeroInflation) 0.02 else 0.0001, b=numeric(p), tmpL=0.25+numeric(sum(m)),
    tmpL1=0.0001+numeric(numb_cor_params), logalpha=2.0, loggammashape=0, u=rep(0,sum(m*q)))


  file_name <- "glmmadmb"

  dat_write(file_name, dat_list)
  pin_write(file_name, pin_list)
  std_file <- paste(file_name, ".std", sep="")
  if(file.exists(std_file))
    file.remove(std_file)

  nbits <- 8 * .Machine$sizeof.pointer
  platform <- if (.Platform$OS.type == "windows") "windows" else {
    ## MacOS detection OK?
    ## http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8497.html
    if (substr(R.version$os,1,6)=="darwin") "macos" else {
      if (R.version$os=="linux-gnu") "linux" else {
        stop("glmmadmb binary not available for this OS")
        ## FIXME: allow user to supply their own binary?
      }
    }
  }
  execname <- if (platform=="windows") paste(file_name,"exe",sep=".") else file_name
  bin_loc <- system.file("bin",paste(platform,nbits,sep=""),
                         execname,
                         package="glmmADMB")
  ## FIXME: for what platforms do we really need to copy the binary?
  ##  can't we just run it in place?  Or does it do something silly and produce
  ##  output in the directory in which the binary lives rather than the
  ##  current working directory (feel like I struggled with this earlier
  ##  but have now forgotten -- ADMB mailing list archives??)
  if (nchar(bin_loc)==0) stop("glmmadmb binary should be available, but isn't")
  if (platform=="windows") {
    cmd <- paste("\"",bin_loc, "\"", " ", cmdoptions, sep="")
    shell(cmd, invisible=TRUE)
  } else {
    file.copy(bin_loc,".")
    Sys.chmod(file_name,mode="0755") ## file.copy strips executable permissions????
    cmd2 <- paste("./", file_name, " ", cmdoptions, sep="")
    sys.result <- system(cmd2,intern=!verbose)
    unlink(file_name)
  }

  if(!file.exists(std_file))
    stop("The function maximizer failed")
  tmp <- read.table(paste(file_name,"std",sep="."), skip=1,as.is=TRUE)
  ## FIXME: could we change the TPL file to write everything out in full precision??
  ##   ... otherwise to read .par file or binary versions ...
  ## BMB: if 'std dev' were written without a space we could use header=TRUE
  tmpindex <- tmp[,2]
  out <- list(n=n, q=q, formula=formula, call=call,
              family=family, corStruct=corStruct, impSamp=impSamp,
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
  } else if (family=="gamma") {
  }
  
  if(!missing(link))
    out$link <- link
  if(has_rand) {
    ## BMB: fixme! make sure this works for multiple random effects
    Svec <- tmp[tmpindex=="S",3]
    parnames <- lapply(REmat$mmat,colnames)
    groupnames <- names(REmat$mmats)
    dn <- mapply(list,
                 parnames,
                 parnames,SIMPLIFY=FALSE)
    out$S <- mapply(matrix,split(Svec,rep(1:length(m),m)),
                    nrow=m,
                    ncol=m,
                    dimnames=dn,
                    MoreArgs=list(byrow=FALSE),
                    SIMPLIFY=FALSE)
    sd_Svec <- tmp[tmpindex=="S",4]
    out$sd_S <- mapply(matrix,split(sd_Svec,rep(1:length(m),m)),
                       nrow=m,
                       ncol=m,
                       dimnames=dn,
                       MoreArgs=list(byrow=FALSE),
                       SIMPLIFY=FALSE)
    names(out$S) <- names(out$sd_S) <- groupnames
    if(corStruct == "diag")
      {
        replace_offdiag <- function(m,r=0) {
          m[upper.tri(m) | lower.tri(m)] <- r
          m
        }
        out$S <- lapply(out$S,replace_offdiag)
        out$sd_S <- lapply(out$sd_S,replace_offdiag)
        ## FIXME: replace with NA for sd_S?
      }
    ## out$random <- random
    ## FIXME: check dimnames for a wider range of cases ...
    uvec <- tmp[tmpindex=="u",3]
    ulist <- split(uvec,rep(1:length(q),q))
    dn <- mapply(list,
                 lapply(REmat$mmat,attr,which="levels"),
                 parnames,
                 SIMPLIFY=FALSE)
    out$U <- mapply(matrix,
                    ulist,
                    ncol=m,
                    nrow=q,
                    dimnames=dn,
                    SIMPLIFY=FALSE)
    names(out$U) <- groupnames
    ## BMB/FIXME: should this be byrow or not? check ...
    ii <- 1
    allU <- matrix(nrow=n,ncol=sum(m))
    for (i in 1:length(m)) {
      allU[,ii:(ii+m[i]-1)] <- out$U[[i]][c(REmat$codes[[i]]),]
      ii <- ii+m[i]
    }
    ## U <- matrix(as.numeric(tmp[tmpindex=="u",3]), ncol=m, byrow=TRUE,
    ## ## dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    ## dimnames=list(NULL,colnames(Z)))
    sd_uvec <- tmp[tmpindex=="u",4]
    sd_ulist <- split(sd_uvec,rep(1:length(q),q))
    out$sd_U <- mapply(matrix,
                       sd_ulist,
                       ncol=m,
                       nrow=q,
                       dimnames=dn,
                       SIMPLIFY=FALSE)
    ## FIXME: check byrow = TRUE or FALSE
    ## matrix(as.numeric(tmp[tmpindex=="u",4]), ncol=m, byrow=TRUE,
    ## dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    ##    dimnames=list(NULL,colnames(Z)))
  } else {
    out$U <- out$sd_U <- matrix(rep(0,q), ncol=1, byrow=TRUE)
  } ## !has_rand

  mu <- as.numeric(X %*% out$b)
  ## BMB: doesn't include influence of random effects?
  lambda <- 0

  ## for(i in 1:n)
  ## lambda[i] <- exp(mu[i] + rowSums(Z * allU))
  if (has_rand) {
    eta <- mu + rowSums(Z*allU)
  } else eta <- mu
  
  lambda <- ilinkfun(eta)

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
  out$npar <- tmp$npar   ## BMB: should this be total number of parameters or number of fixed parameters?
  bpar <- bar_read(file_name,n=tmp$npar+1)[-1] ## ZI parameter comes first
  tmp$beta <- bpar
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

  ## BMB: warning/convergence code for bad hessian?

  ## FIXME: figure out names of parameters
  if (mcmc) {
    out$mcmc <- read_psv(file_name)
  }
  
  class(out) <- "glmm.admb"

  return(out)
}
