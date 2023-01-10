predict_ff_regression <- function(object, Xnew){

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

  gpY <- object$model.details$gpY
  Bhat <- object$model.details$Bhat
  np_test <- length(Xnew)
  n_test <- dim(Xnew[[1]])[1]
  p <- dim(Y)[2]

  comp_Y <- object$fpca.results$Y$PCAcoef
  mean_Y <- object$fpca.results$Y$meanScore
  obj <- object$fpca.results$X

  sco_X_test <- list()
  for(i in 1:np_test)
    sco_X_test[[i]] <- getPCA.test(object = obj[[i]], data = Xnew[[i]])

  Xkf <- do.call(cbind, sco_X_test)

  Ypreds <- matrix(, nrow = n_test, ncol = p)
  for(k in 1:n_test){
    Xk = Xkf[k,]
    model_k <- pred.fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
    Ypreds[k,] <- eval.fd(gpY, model_k)
  }

  return(Ypreds)
}
