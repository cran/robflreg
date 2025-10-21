pqr <- function(y, x, tau, h, qc.type){
  x0 <- x
  y0 <- y
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(y)[2]
  
  x <- scale(x, scale = F)
  y <- scale(y, scale = F)
  
  xw <- matrix(, p, h)
  yw <- matrix(, q, h)
  xs <- matrix(, n, h)
  ys <- matrix(, n, h)
  xl <- matrix(, p, h)
  yl <- matrix(, q, h)
  
  for(i in 1:h)
  {
    qcv <- matrix(, p, q)
    for(qi in 1:q)
      qcv[,qi] <- qcov(y[,qi], x, tau, qc.type)
    qsvd <- svd(qcv)
    
    xw. <- qsvd$u[,1]
    yw. <- qsvd$v[1,]
    
    buff <- svd_flip(xw., yw.)
    xw. <- buff$x
    yw. <- buff$y
    
    xs. <- x %*% xw.
    ys. <- y %*% yw.
    
    xl. <- c(t(x) %*% xs.) / c(t(xs.) %*% xs.)
    x <- x - xs. %*% xl.
    yl. <- c(t(xs.) %*% y) / c(t(xs.) %*% xs.)
    y <- y - xs. %*% yl.
    
    xw[,i] <- xw.
    yw[,i] <- yw.
    xs[,i] <- xs.
    ys[,i] <- ys.
    xl[,i] <- xl.
    yl[,i] <- yl.
  }
  
  xr <- xw %*% ginv(t(xl) %*% xw)
  yr <- yw %*% ginv(t(yl) %*% yw)
  
  fin.cf <- matrix(, (p+1), q)
  pqr.cf <- matrix(, (h+1), q)
  for(ic in 1:q){
    model.ic <- rq(y0[,ic]~xs, tau)
    temp <- xr %*% as.matrix(model.ic$coefficients[-1])
    pqr.cf[,ic] <- model.ic$coefficients
    fin.cf[,ic] <- c(model.ic$coefficients[1], temp)
  }
  
  return(list(T = xs, R = xr, P = xl, W = xw, d.coef = fin.cf,
              pqr.coef = pqr.cf))
}