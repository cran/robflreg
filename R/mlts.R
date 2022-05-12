mlts <-
function(X, Y, gamma, ns = 500, nc = 10, delta = 0.01){
  
  d <- dim(X)
  n <- d[1]
  p <- d[2]
  q <- ncol(Y)
  h <- floor(n*(1-gamma))+1
  obj0 <- 1e10
  for(i in 1:ns){
    sorted <- sort(runif(n), na.last = NA, index.return=TRUE)
    istart <- sorted$ix[1:(p+q)]
    xstart <- X[istart,]
    ystart <- Y[istart,]
    bstart <- solve(t(xstart) %*% xstart,t(xstart) %*% ystart)
    sigmastart <- (t(ystart-xstart %*% bstart)) %*% (ystart-xstart%*%bstart)/q
    for(j in 1:nc){
      res <- Y - X %*% bstart
      tres <- t(res)
      dist2 <- colMeans(solve(sigmastart,tres)*tres)
      sdist2 <- sort(dist2, na.last = NA, index.return = TRUE)
      idist2 <- sdist2$ix[1:h]
      xstart <- X[idist2,]
      ystart <- Y[idist2,]
      bstart <- solve(t(xstart) %*% xstart,t(xstart) %*% ystart)
      sigmastart <- (t(ystart-xstart %*% bstart)) %*% (ystart-xstart %*% bstart)/(h-p)
    }
    obj <- det(sigmastart)
    if(obj < obj0){
      result.beta <- bstart
      result.sigma <- sigmastart
      obj0 <- obj
    }
  }
  cgamma <- (1-gamma)/pchisq(qchisq(1-gamma,q),q+2)
  result.sigma <- cgamma * result.sigma
  res <- Y - X %*% result.beta
  tres <- t(res)
  result.dres <- colSums(solve(result.sigma,tres)*tres)
  result.dres <- sqrt(result.dres)
  
  qdelta <- sqrt(qchisq(1-delta,q))
  good <- (result.dres <= qdelta)
  xgood <- X[good,]
  ygood <- Y[good,]
  result.betaR <- solve(t(xgood) %*% xgood,t(xgood) %*% ygood)
  
  list(betaR=t(result.betaR))
}
