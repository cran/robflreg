get.ff.coeffs <-
function(object){
  
  vars <- object$model.details$var.used
  np <- length(vars)
  ncX <- object$model.details$ncompX
  gpY <- object$model.details$gpY
  gpX <- object$model.details$gpX
  Bhat <- object$model.details$Bhat
  evaly <- object$fpca.results$Y$evalbase
  fcomp_Y <- object$fpca.results$Y$PCAcoef$coefs
  evalx <- list()
  fcomp_X <- list()
  for(i in 1:np){
    evalx[[i]] <- object$fpca.results$X[[i]]$evalbase
    fcomp_X[[i]] <- object$fpca.results$X[[i]]$PCAcoef$coefs
  }
  
  coeffs = list()
  km = 1
  for(j in 1:np){
    coeffs[[j]] = evalx[[j]] %*% (fcomp_X[[j]] %*% Bhat[km: (km+ncX[j]-1),] %*% t(fcomp_Y)) %*% t(evaly)
    km = j*ncX[j]+1
  }
  
  return(list(vars = vars,
              gpY = gpY,
              gpX = gpX,
              coefficients = coeffs))
    

  
}
