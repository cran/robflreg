fpqr <- function(y, x, h=NULL, tau, nby=NULL, nbx=NULL, gpy=NULL, gpx=NULL, 
                 qc.type = c("dodge","choi","li"), hs=1:5, nbys = c(4,8,10,20),
                 nbxs = c(4,8,10,20), nfold=5, CV = TRUE){
  
  if(!is.matrix(x))
    stop("Error!! x must be a matrix!")
  if(!is.matrix(y))
    stop("Error!! y must be a matrix!")
  if(!is.null(gpy)){
    if(length(gpy) != dim(y)[2])
      stop("Error!! The number of columns of y must be equal to the length of grid points: gpy!")
  }
  if(!is.null(gpx)){
    if(length(gpx) != dim(x)[2])
      stop("Error!! The number of columns of x must be equal to the length of grid points: gpx!")
  }
  
  if(is.null(gpx))
    gpx <- seq(0, 1, length.out = dim(x)[2])
  if(is.null(gpy))
    gpy <- seq(0, 1, length.out = dim(x)[2])
  
  qc.type <- match.arg(qc.type)
  
  if(CV == TRUE){
    optimod <- paramCV(y=y, x=x, hs = hs, tau = tau, nbys = nbys, nbxs = nbxs, 
                       gpy = gpy, gpx = gpx, nfold = nfold, qc.type = qc.type)
    h <- optimod$h
    nby <- optimod$nby
    nbx <- optimod$nbx
  }

  BS.sol.x <- getAmat(data = x, nbf = nbx, gp = gpx)
  x1 <- BS.sol.x$Amat
  BS.sol.y <- getAmat(data = y, nbf = nby, gp = gpy)
  y1 <- BS.sol.y$Amat

  m.pqr <- pqr(y = y1, x = x1, tau = tau, h = h, qc.type = qc.type)
  fits <- (cbind(1,x1) %*% m.pqr$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)

  b.hat <- BS.sol.x$evalbase %*% solve(BS.sol.x$sinp_mat) %*% m.pqr$d.coef[-1,] %*%
    solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  b0.hat <- t(solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)) %*% as.matrix(m.pqr$d.coef[1,])

  model.details <- list()
  model.details$nbx <- nbx
  model.details$gpx <- gpx
  model.details$m.pqr <- m.pqr
  model.details$BS.sol.x <- BS.sol.x
  model.details$BS.sol.y <- BS.sol.y

  return(list(fitted.values = fits, b0.hat = b0.hat, b.hat = b.hat, mdts = model.details))

}
