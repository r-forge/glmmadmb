"glmm.admb" <- function(fixed, random, group, data, family="poisson", link, corStruct="diag", impSamp=0, easyFlag=TRUE,
                        zeroInflation=FALSE, imaxfn=10, save.dir=NULL, offset)
{
  olddir <- getwd()
  dirname <- ifelse(is.null(save.dir), "_glmm_ADMB_temp_dir_", save.dir)
  dir.create(dirname)
  setwd(dirname)

  if(is.null(save.dir))
    on.exit(unlink(dirname, recursive=TRUE))
  call <- match.call()
  if(!(corStruct == "diag" | corStruct == "full"))
    stop("Argument \"corStruct\" must be either \"diag\" or \"full\"")
  if(!(family == "poisson" | family == "nbinom" | family ==
       "binomial"))
    stop("Argument \"family\" must be either \"poisson\", \"binomial\"  or \"nbinom\"")
  if(missing(random) & ((!missing(impSamp))|(!missing(corStruct))))
    stop("When \"random\" is missing neither of \"impSamp\" or \"corStruct\" make sense")
  if(missing(group) || (!(is.character(group)&(length(group)==1))) || !member(group,names(data)))
    stop("Argument \"group\" must be a character string spesifying the name of the grouping variable (also when \"random\" is missing)")
  if(missing(link) & family=="binomial")
    stop("Argument \"link\" must be provided for the \"binomial\" family)")
  if((!missing(offset)) && ((!is.character(offset))||(length(offset)!=1)||!member(offset,names(data))))
    stop("An error occured in relatation to the argument \"offset\". It must be a character string spesifying of the variable holding the offset")
  if(!missing(offset))
    Offset = data[[offset]]
  else
    Offset = rep(0,nrow(data))

  tmpu <- table(data[[group]])
  tmpu[] <- 1:length(tmpu)
  tmpI <- tmpu[as.character(data[[group]])]
  data <- data[order(tmpI), ]
  n <- nrow(data)
  y <- data[[as.character(fixed[2])]]
  X <- model.matrix(fixed, data)
  p <- ncol(X)
  II <- as.numeric(tmpu[as.character(data[[group]])])
  group_d <- c(1, (2:n)[diff(II)==1], n+1)
  if(missing(random))
    Z <- model.matrix(~1, data)
  else
    Z <- model.matrix(random, data)
  m <- ncol(Z)
  q <- length(unique(II))

  cmdoptions = paste("-maxfn 500", ifelse(impSamp==0,"",paste("-is",impSamp)))
  dat_list = list(n=n, y=y, p=p, X=X, q=q, m=m, Z=Z, group_d=group_d, II=II, cor_flag=ifelse(corStruct=="full",1,0),
    like_type_flag=as.numeric(family=="poisson"||((family == "binomial")&&(link=="logit"))),
    no_rand_flag=as.numeric(missing(random)), easy_flag=as.numeric(easyFlag), zi_flag=as.numeric(zeroInflation),
    intermediate_maxfn=10, offset=Offset)
  pin_list=list(pz=ifelse(zeroInflation,0.02,1e-04), b=numeric(p), tmpL=0.25+numeric(m), tmpL1=1e-04+numeric(m*(m-1)/2),
    logalpha=2.0, kkludge=0, u=matrix(numeric(q*m),q,m))

  if(family == "binomial")
  {
    if(!all(member(y, 0:1)))
      stop("Only Bernoulli responce (0 or 1) allowed currently")
    file_name = "bvprobit"
    dat_list = dat_list[c("n", "y", "p", "X", "q", "m", "Z", "group_d", "cor_flag", "like_type_flag", "no_rand_flag")]
    pin_list = pin_list[c("b", "tmpL", "tmpL1", "kkludge", "u")]
  }
  else
  {
    file_name = "nbmm"
    dat_list = dat_list[c("n", "y", "p", "X", "q", "m", "Z", "II", "cor_flag", "easy_flag", "zi_flag", "like_type_flag",
      "no_rand_flag", "intermediate_maxfn", "offset")]
    pin_list = pin_list[c("pz", "b", "tmpL", "tmpL1", "logalpha", "kkludge", "u")]
    if(!missing(link))
      warning("Argument \"link\" ignored for other familes than \"binomial\" (exponential link used)")
  }

  dat_write(file_name, dat_list)
  pin_write(file_name, pin_list)
  std_file = paste(file_name, ".std", sep="")
  file.remove(std_file)

  if(.Platform$OS.type == "windows")
  {
    shell(paste(.path.package("glmmADMB"),"/admb/",file_name,".exe"," ",cmdoptions,sep=""), invisible=TRUE)
  }
  else
  {
    system(paste("cp ",.path.package("glmmADMB"),"/admb/",file_name," .",sep=""))
    system(paste("./",file_name," ",cmdoptions,sep=""))
    unlink(file_name)
  }

  if(!file.exists(std_file))
    stop("The function maximizer failed")
  tmp <- read.table(paste(file_name,".std",sep=""), skip=1)
  tmpindex <- as.character(tmp[, 2])
  out <- list(n=n, q=q, fixed=fixed, call=call, group=group, family=family, corStruct=corStruct, impSamp=impSamp,
              easyFlag=easyFlag, zeroInflation=zeroInflation)
  if(zeroInflation)
    out$pz=as.numeric(tmp[tmpindex=="pz", 3])
  out$b = as.numeric(tmp[tmpindex=="real_b", 3])
  names(out$b) = colnames(X)
  out$stdbeta = as.numeric(tmp[tmpindex=="real_b", 4])
  names(out$stdbeta) = colnames(X)

  if(family == "nbinom")
  {
    out$alpha <- as.numeric(tmp[tmpindex=="alpha", 3])
    out$sd_alpha <- as.numeric(tmp[tmpindex=="alpha", 3])
  }
  if(!missing(random))
  {
    out$S = matrix(as.numeric(tmp[tmpindex=="S",3]), nrow=m, dimnames=list(colnames(Z),colnames(Z)))
    out$sd_S = matrix(as.numeric(tmp[tmpindex=="S",4]), nrow=m, dimnames=list(colnames(Z), colnames(Z)))
    if(corStruct == "diag")
    {
      out$S[below(m,TRUE)|t(below(m,TRUE))] = 0
      out$sd_S[below(m,TRUE)|t(below(m,TRUE))] = 0
    }
    out$random = random
  }
  if(!missing(link))
    out$link = link
  if(!missing(random))
  {
    U = matrix(as.numeric(tmp[tmpindex=="u",3]), ncol=m, byrow=TRUE,
      dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
    out$U = U
    out$sd_U = matrix(as.numeric(tmp[tmpindex=="u",4]), ncol=m, byrow=TRUE,
      dimnames=list(data[,group][group_d[-1]-1],colnames(Z)))
  }
  else
    U = matrix(rep(0,q), ncol=1, byrow=TRUE)

  mu <- as.numeric(X %*% out$b)
  lambda <- 0
  for(i in 1:n)
    lambda[i] <- exp(mu[i] + sum(Z[i, ]*U[II[i],]))
  if(family == "binomial")
    out$fitted <- lambda / (1+lambda)
  else
    out$fitted <- lambda

  tmpsd <- switch(family,
                  poisson=sqrt(lambda),
                  nbinom=sqrt(lambda*(1+lambda/out$alpha)),
                  binomial=sqrt(out$fitted*(1-out$fitted)))
  out$residuals <- as.numeric(y-lambda) / tmpsd
  tmp = par_read(file_name)
  out$npar <- as.numeric(scan(paste(file_name,".par",sep=""), what="", quiet=TRUE)[6])
  out$loglik = tmp$loglik
  out$gradloglik = tmp$gradient

  if(abs(out$gradloglik) >= 0.001)
    warning("Proper convergence could not be reached")
  setwd(olddir)
  class(out) <- "glmm.admb"

  return(out)
}
