glmmadmb <- function(fixed, random, data, family="poisson", link, impSamp=0, easyFlag=TRUE,
                      zeroInflation=FALSE, imaxfn=10, save.dir=NULL, offset, verbose=FALSE)
{
#  dirname = if(is.null(save.dir)) "_glmm_ADMB_temp_dir_" else save.dir
#  if(!file_test("-d",dirname))
#    dir.create(dirname)
#  owd = setwd(dirname); on.exit(setwd(owd))
#  if(is.null(save.dir))
#    on.exit(unlink(dirname,recursive=TRUE), add=TRUE)

  # Checking command line arguments
  call = match.call()
  if(!(family %in% c("poisson","nbinom","binomial")))
    stop("Argument \"family\" must be either \"poisson\", \"binomial\"  or \"nbinom\"")
  if(missing(random) && (!missing(impSamp) || !missing(corStruct)))
    stop("When \"random\" is missing, neither of \"impSamp\" or \"corStruct\" make sense")
  if(missing(link) && family=="binomial")
    stop("Argument \"link\" must be provided for the \"binomial\" family)")
  if(!missing(offset) && (!is.character(offset) || length(offset)!=1 || !(offset %in% names(data))))
    stop("\"offset\" must be a character string specifying the variable holding the offset")

  # Offset handling (may need revision?)
  has_offset = !missing(offset)
  Offset = if(has_offset) data[[offset]] else rep(0,nrow(data))

  # Extract data
  y = data[[as.character(fixed[2])]]
  n = length(y)
  X = model.matrix(fixed, data)
  p = ncol(X)

  REmat = process_randformula(random,data=data)

  m = as.numeric(lapply(REmat$mmats, function(x) ncol(x)))		# Number of random effects within each crossed term
  M = length(m)		 						# Number of crossed terms
  q = as.numeric(lapply(REmat$codes, function(x) length(unique(x))))	# Number of levels of each crossed term

  # Construct design matrix (Z) for random effects
  for(i in 1:M)
  {
    tmp = REmat$mmats[[i]]
    colnames(tmp) = paste(names(REmat$mmats)[i],colnames(tmp),sep=":")
    if(i==1)
      Z = tmp
    else
      Z = cbind(Z,tmp)
  }


  # Make matrix of pointers into ADMB vector of random effects "u".
  # Rows in II contains pointers for a given data point 
  # First offsets the collumns of II to point into different segments of u
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

  # Splits u into "correlation blocks" (yet no way of specifying uncorrelated random effects
  cor_block_start = cumsum(m)-m[1]+1
  cor_block_stop = cumsum(m)
  numb_cor_params = max(sum(m*(m-1)/2),1)	# For tecnical reason we need at least 1 parameter: see glmmadmb.tpl
  
  cmdoptions = paste("-maxfn 500", if(impSamp==0) "" else paste("-is",impSamp))
  dat_list = list(n=n, y=y, p=p, X=X, M=M, q=q, m=m, ncolZ=ncol(Z),Z=Z, II=II, 
		   	cor_flag=rep(0,M),
			cor_block_start=cor_block_start,
		   	cor_block_stop=cor_block_stop,
			numb_cor_params=numb_cor_params,
                   	like_type_flag=as.numeric(family=="poisson"||(family=="binomial"&&link=="logit")),
                   	no_rand_flag=as.numeric(missing(random)), 
                   	zi_flag=as.numeric(zeroInflation), 
			intermediate_maxfn=10, 
			has_offset=as.numeric(has_offset), 
			offset=Offset)
  pin_list = list(pz=if(zeroInflation) 0.02 else 0.0001, b=numeric(p), tmpL=0.25+numeric(sum(m)),
                   tmpL1=0.0001+numeric(numb_cor_params), logalpha=2.0, u=rep(0,sum(m*q)))


  file_name = "glmmadmb"

  dat_write(file_name, dat_list)
  pin_write(file_name, pin_list)
  std_file = paste(file_name, ".std", sep="")
  if(file.exists(std_file))
    file.remove(std_file)

  system("./glmmadmb -shess")


}
