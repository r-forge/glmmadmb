".First.lib" <- function(lib, pkg)
{
  cat("\n\n")
  cat(paste(rep("=",options()$width),collapse=""), "\n\n")
  cat("Welcome to glmmADMB\n")
  if(.Platform$OS.type == "windows")
  {
    cat("\n")
    cat("WINDOWS USERS: Please press (Ctrl+W) or use the menu \n")
    cat("(Misc->Buffered output) to switch the delayed output off.\n")
  }
  cat("\n")
  cat(paste(rep("=",options()$width),collapse=""), "\n\n")
}
