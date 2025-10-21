predict_fpqr <- function(object, xnew){
  
  nbx <- object$mdts$nbx
  gpx <- object$mdts$gpx
  m.pqr <- object$mdts$m.pqr
  BS.sol.x <- object$mdts$BS.sol.x
  BS.sol.y <- object$mdts$BS.sol.y
  
  Bs.sol.xnew <- getAmat(data = xnew, nbf = nbx, gp = gpx)
  x1new <- Bs.sol.xnew$Amat
  
  preds <- (cbind(1,x1new) %*% m.pqr$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
  
  return(preds)
  
}