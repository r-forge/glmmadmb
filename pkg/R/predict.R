nobs.glmmadmb <- function(object,...) {
  nrow(object$frame)
}

predict.glmmadmb <- function(object, newdata=NULL, 
                    type=c("link","response"), se.fit=FALSE, ...) {
  if (se.fit && type=="response") {
    warning("se.fit && type='response': setting se.fit to NA")
  }
  type <- match.arg(type)
  ## Construct model matrix, nobs x np
  if (is.null(newdata)) newdata <- object$frame
  form <- as.formula(as.character(object$fixed)[-2])
  X <- model.matrix(form, data=newdata)
  beta <- as.vector(object$b)
  phat <- X %*% beta
  offset <- rep(0, nrow(X))
  tt <- object$terms
  ## copied from predict.lm, unpacked slightly for ease of debugging
  off.num <- attr(tt, "offset")
  if (!is.null(off.num))  {
    for (i in off.num) {
      cur.offset <- eval(attr(tt, "variables")[[i]], newdata)
      offset <- offset + cur.offset
    }
  }
  phat <- c(phat + offset)
  if (se.fit) {
    stderr <- c(sqrt(diag(X %*% vcov(X) %*% t(X))))
  }
  if (type=="response") {
    phat <- object$ilinkfun(phat)
    stderr <- NA
  }
  if (se.fit) list(fit=phat,se.fit=stderr) else phat
}
