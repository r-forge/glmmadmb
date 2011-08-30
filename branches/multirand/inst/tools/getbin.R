RELEASE <- "r135"
w7ver <- "mingw"
buildbot_base <- "http://www.admb-project.org/buildbot/glmmadmb/glmmadmb-admb-"
## assume are starting from the root of the package directory
platform_str <- c(linux32="ubuntu11-gcc4.5-32bit-",
                  linux64="ubuntu11-gcc4.5-64bit-",
                  macos32="macos10.6-xcode3.2-32bit-",
                  macos64="macos10.6-xcode3.2-64bit-",
                  windows32="windows7-borland-32bit-",
                  windows64="windows7-vc10-64bit-")
  ## c(linux32="linux-gcc4.5.2-32bit",
  ##   linux64="linux-gcc4.5.2-64bit",
  ##   macos32="macos10.6.7-xcode3.2.6-32bit",
  ##   macos64="macos10.6.7-xcode3.2.6-64bit")
if (!file.exists("inst/bin")) {
  setwd("../..")  ## try moving up
  if (!file.exists("inst/bin")) {
    stop("can't find bin directory -- are you in the package root?")
  }
}
setwd("inst/bin")
for (i in names(platform_str)) {
  setwd(i)
  srcext <- if (grepl("windows",i)) ".exe" else ".bin"
  destext <- if (grepl("windows",i)) ".exe" else ""
  fn <- paste(buildbot_base,platform_str[i],RELEASE,
                      srcext,sep="")
  download.file(fn,
                destfile=paste("glmmadmb",destext,sep=""))
  setwd("..")
}
  
## glmmadmb-r107-bcc5.5-32bit 	0B 	[text/html] 	
## glmmadmb-r107-linux-gcc4.4.3-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-64bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-solarisstudio12-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-64bit 	2M 	[text/html] 
