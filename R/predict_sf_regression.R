predict_sf_regression <-
function(object, Xnew){
  
  if(!is.list(Xnew))
    stop("Error!! Xnew must be a list!")
  if(sum(sapply(1:length(Xnew),function(i){!is.matrix(Xnew[[i]])})))
    stop("Error!! Xnew must be a list and all of its components must be matrix!")
  dimXnew <- sapply(1:length(Xnew),function(i){dim(Xnew[[i]])[1]})
  if((length(unique(dimXnew))!=1))
    stop("Error!! All the components of Xnew must be matrix and they must have the same number of rows!")
  
  Y <- object$data$Y
  X <- object$data$X
  
  dimXnew2 <- sapply(1:length(Xnew),function(i){dim(Xnew[[i]])[2]})
  dimX2 <- sapply(1:length(X),function(i){dim(X[[i]])[2]})
  if(!all(dimXnew2 == dimX2))
    stop("Error!! All opposing elements of X and Xnew must have the same number of columns!")
  if(length(X) != length(Xnew))
    stop("Error!! Xnew and X must have the same number of variables!")
  
  nbasis <- object$model.details$nbasis
  gp <- object$model.details$gp
  ncomp <- object$model.details$ncomp
  emodel <- object$model.details$emodel
  Bhat <- object$model.details$Bhat
  np <- length(X)
  np_test <- length(Xnew)
  n_test <- dim(Xnew[[1]])[1]

  mean_Y <- object$model.details$mean_Y
  
  sco_X_test <- list()
  for(i in 1:np_test)
    sco_X_test[[i]] <- getPCA.test(data = Xnew[[i]], bs_basis = object$fpca.results[[i]]$bs_basis,
                                   PCAcoef = object$fpca.results[[i]]$PCAcoef,
                                   gp = gp[[i]], emodel = emodel)
  
  Ypreds <- do.call(cbind, sco_X_test) %*% Bhat + mean_Y
  
  return(Ypreds)
}
