predict_sf_regression <-
function(object, Xnew, Xnew.scl = NULL){

  if(!is.list(Xnew))
    stop("Error!! Xnew must be a list!")
  if(sum(sapply(1:length(Xnew),function(i){!is.matrix(Xnew[[i]])})))
    stop("Error!! Xnew must be a list and all of its components must be matrix!")
  dimXnew <- sapply(1:length(Xnew),function(i){dim(Xnew[[i]])[1]})
  if((length(unique(dimXnew))!=1))
    stop("Error!! All the components of Xnew must be matrix and they must have the same number of rows!")
  if(!is.null(Xnew.scl)){
    if(!is.matrix(Xnew.scl))
      Xnew.scl <- as.matrix(Xnew.scl)
  }
  if(!is.null(Xnew.scl)){
    if(dim(Xnew.scl)[1] != dim(Xnew[[1]])[1])
      stop("Error!! Functional and predictor variables must have the same number of rows!")
  }

  Y <- object$data$Y
  X <- object$data$X
  X.scl <- object$data$X.scl

  dimXnew2 <- sapply(1:length(Xnew),function(i){dim(Xnew[[i]])[2]})
  dimX2 <- sapply(1:length(X),function(i){dim(X[[i]])[2]})
  if(!all(dimXnew2 == dimX2))
    stop("Error!! All opposing elements of X and Xnew must have the same number of columns!")
  if(length(X) != length(Xnew))
    stop("Error!! Xnew and X must have the same number of variables!")
  if(!is.null(Xnew.scl)){
    if(dim(Xnew.scl)[2] != dim(X.scl)[2])
      stop("Error!! Xnew.scl and X.scl must have the same number of variables!")
  }

  emodel <- object$model.details$emodel
  Bhat <- object$model.details$Bhat
  Bshat <- object$model.details$Bshat
  np_test <- length(Xnew)

  mean_Y <- object$model.details$mean_Y
  obj <- object$fpca.results

  sco_X_test <- list()
  for(i in 1:np_test)
    sco_X_test[[i]] <- getPCA.test(object = obj[[i]], data = Xnew[[i]])

  if(!is.null(Xnew.scl)){
    if(emodel == "classical")
      mean_Xnew.scl <- apply(Xnew.scl, 2, mean)
    if(emodel == "robust")
      mean_Xnew.scl <- apply(Xnew.scl, 2, median)
    Xnew.scl <- scale(Xnew.scl, center = mean_Xnew.scl, scale = F)
  }

  if(!is.null(Xnew.scl)){
    Xall <- cbind(Xnew.scl, do.call(cbind, sco_X_test))
    Ball <- c(Bshat, Bhat)
  }else{
    Xall <- do.call(cbind, sco_X_test)
    Ball <- Bhat
  }

  Ypreds <- Xall %*% Ball + mean_Y

  return(Ypreds)
}
