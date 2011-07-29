RELEASE <- "r107"
buildbot_base <- "http://www.admb-project.org/buildbot/glmmadmb/glmmadmb-"
## assume are starting from the root of the package directory
platform_str <- c(linux32="linux-gcc4.5.2-32bit",
                  linux64="linux-gcc4.5.2-64bit",
                  macos32="macos10.6.7-xcode3.2.6-32bit",
                  macos64="macos10.6.7-xcode3.2.6-64bit")
if (!file.exists("inst/bin"))
  stop("can't find bin directory -- are you in the package root?")
setwd("inst/bin")
for (i in names(platform_str)) {
  setwd(i)
  download.file(paste(buildbot_base,RELEASE,"-",platform_str[i],sep=""),
                destfile="glmmadmb")
  setwd("..")
}
  
## glmmadmb-r107-bcc5.5-32bit 	0B 	[text/html] 	
## glmmadmb-r107-linux-gcc4.4.3-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-32bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-gcc4.5.2-64bit 	1M 	[text/html] 	
## glmmadmb-r107-linux-solarisstudio12-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-32bit 	2M 	[text/html] 	
## glmmadmb-r107-macos10.6.7-xcode3.2.6-64bit 	2M 	[text/html] 
