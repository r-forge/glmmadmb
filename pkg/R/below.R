"below" <- function(n, strictly=FALSE)
{
  M <- matrix(TRUE, n, n)

  M[rep(1:n,n)<rep(1:n,rep(n,n))] <- FALSE
  if(strictly)
    diag(M) = FALSE

  M
}
