mcd <-
function(X, Y){
  
  p <- dim(X)[2]
  q <- dim(Y)[2]

  z <- cbind(Y, X)
  z2 <- covMcd(z, raw.only = FALSE, alpha = 0.75)
  At <- z2$cov

  a <- (At[1:q, (q+1):(p+q)]) %*% ginv(At[(q+1):(p+q), (q+1):(p+q)])
  
  
  return(list(betaR = a))
}
