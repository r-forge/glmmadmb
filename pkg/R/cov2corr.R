"cov2corr" <-
function (m)
diag(1/sqrt(diag(m))) %*% m %*% diag(1/sqrt(diag(m)))
