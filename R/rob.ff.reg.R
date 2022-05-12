rob.ff.reg <-
function(Y, X, model = c("full", "selected"), emodel = c("classical", "robust"),
                      fmodel = c("MCD", "MLTS", "MM", "S", "tau"), nbasisY = NULL, nbasisX = NULL, 
                      gpY = NULL, gpX = NULL, ncompY = NULL, ncompX = NULL){
  
  if(!is.matrix(Y))
    stop("Error!! Y must be matrix!")
  if(!is.list(X))
    stop("Error!! X must be a list!")
  if(sum(sapply(1:length(X),function(i){!is.matrix(X[[i]])})))
    stop("Error!! X must be a list and all of its components must be matrix!")
  dimX <- sapply(1:length(X),function(i){dim(X[[i]])[1]})
  if((length(unique(dimX))!=1))
    stop("Error!! All the components of X must be matrix and they must have the same number of rows!")
  if(dimX[1] != dim(Y)[1])
    stop("Error!! Response and predictor variables must have the same number of rows!")
  if(!is.null(gpY)){
    if(length(gpY) != dim(Y)[2])
      stop("Error!! The number of columns of Y must be equal to the length of grid points (gpY)!")
  }
  if(!is.null(gpX))
  {
    if(!is.list(gpX))
      stop("Error!! gpX must be a list!")
    lengpX <- sapply(1:length(gpX), function(i){length(gpX[[i]])})
    ncX <- sapply(1:length(X), function(i){dim(X[[i]])[2]})
    if(length(gpX) != length(X))
      stop("Error!! The lengths of X and gpX must be equal!")
    if(!all(lengpX == ncX))
      stop("Error!! The number of columns of each component of X must be equal to the length of corresponding component of gpX!")
  }
  if(!is.null(nbasisY)){
    if(!nbasisY >= 4)
      stop("Error!! nbasisY must be greater than three to apply cubic B-spline basis expansion!")
  }
  if(!is.null(nbasisX)){
    if(!is.numeric(nbasisX))
      stop("Error!! nbasisX must be a numeric vector!")
    if(!all(nbasisX >= 4))
      stop("Error!! Each component of nbasisX must be greater than three to apply cubic B-spline basis expansion!")
    if(length(nbasisX) != length(X))
      stop("Error!! The lengths of X and nbasisX mut be equal!")
  }
  if(!is.null(ncompY)){
    if(!ncompY > 0)
      stop("Error!! The number of principal components (ncomp) must be greater than or equal to one!")
  }
  if(!is.null(ncompX)){
    if(!is.numeric(ncompX))
      stop("Error!! ncompX must be a numeric vector!")
    if(!all(ncompX > 0))
      stop("Error!! Each component of ncompX (the numbers of principal components) must be greater than or equal to one!")
    if(length(ncompX) != length(X))
      stop("Error!! The lengths of X and ncompX must be equal!")
  }
  if(!model %in% c("full", "selected"))
    stop("Error!! model must be one of the followings: full or selected!")
  if(!emodel %in% c("classical", "robust"))
    stop("Error!! emodel must be one of the followings: classical or robust!")
  if(emodel == "robust")
  {
    if(!fmodel %in% c("MCD", "MLTS", "MM", "S", "tau"))
      stop("Error!! model must be one of the followings: MCD, MLTS, MM, S, and tau!")
  }
  
  model <- match.arg(model)
  emodel <- match.arg(emodel)
  if(emodel != "classical"){
    fmodel <- match.arg(fmodel)
  }else{
    fmodel <- NULL
  }
  
  if(model == "full"){
    
    np <- length(X)
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    px <- numeric()
    for(ip in 1:np)
      px[ip] <- dim(X[[ip]])[2]

    if(is.null(nbasisY)) nbasisY = 20
    if(is.null(nbasisX)) nbasisX = rep(20, np)
    if(is.null(gpY)) gpY = seq(0, 1, length.out = p)
    if(is.null(gpX)){
      gpX <- list()
      for(iip in 1:np)
        gpX[[iip]] = seq(0, 1, length.out = px[iip])
    }
    if(is.null(ncompY)) ncompY = 4
    if(is.null(ncompX)) ncompX = rep(4, np)
    

    PCA_Y <- getPCA(data = Y, nbasis = nbasisY, ncomp = ncompY,
                        gp = gpY, emodel = emodel)
    sco_Y <- PCA_Y$PCAscore
    comp_Y <- PCA_Y$PCAcoef
    mean_Y <- PCA_Y$meanScore
    
    PCA_X <- list()
    sco_X <- list()
    for(fij in 1:np){
      PCA_X[[fij]] <- getPCA(data = X[[fij]], nbasis = nbasisX[fij], ncomp = ncompX[fij],
                          gp = gpX[[fij]], emodel = emodel)
      sco_X[[fij]] <- PCA_X[[fij]]$PCAscore
    }
    
    Bhat <- est.fun(sco_Y = sco_Y, sco_X = do.call(cbind, sco_X),
                    emodel = emodel, fmodel = fmodel)
    
    Yfit <- matrix(NA, nrow = n, ncol = p)
    for(k in 1:n){
      Xk <- do.call(cbind, sco_X)[k,]
      model_k <- pred.fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
      Yfit[k,] <- eval.fd(model_k, seq(gpY[1], gpY[p], length.out = p))
    }

    resids <- Y - Yfit
    var_index <- c(1:np)
  }
  
  if(model == "selected"){
    X1 <- X
    np1 <- length(X)
    if(np1 == 1) stop("Error!! Number of functional predictors should be greater than one to use variable selection procedure!")
    px <- numeric()
    for(ip in 1:np1)
      px[ip] <- dim(X1[[ip]])[2]
    p <- dim(Y)[2]
    
    if(is.null(nbasisY)) nbasisY = 20
    if(is.null(nbasisX)) nbasisX = rep(20, np1)
    if(is.null(gpY)) gpY = seq(0, 1, length.out = p)
    if(is.null(gpX)){
      gpX <- list()
      for(iip in 1:np1)
        gpX[[iip]] = seq(0, 1, length.out = px[iip])
    }
    if(is.null(ncompY)) ncompY = 4
    if(is.null(ncompX)) ncompX = rep(4, np1)
    
    svar <- var.sel(Y = Y, X = X, nbasisY = nbasisY, nbasisX = nbasisX,
                    gpY = gpY, gpX = gpX, emodel = emodel)
    mindex <- svar$maine

    X <- X[mindex]

    np <- length(X)
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    px <- px[mindex]
    nbasisX <- nbasisX[mindex]
    gpX <- gpX[mindex]
    ncompX <- ncompX[mindex]
    
    
    PCA_Y <- getPCA(data = Y, nbasis = nbasisY, ncomp = ncompY,
                        gp = gpY, emodel = emodel)
    sco_Y <- PCA_Y$PCAscore
    comp_Y <- PCA_Y$PCAcoef
    mean_Y <- PCA_Y$meanScore
    
    PCA_X <- list()
    sco_X <- list()
    for(fij in 1:np){
      PCA_X[[fij]] <- getPCA(data = X[[fij]], nbasis = nbasisX[fij], ncomp = ncompX[fij],
                          gp = gpX[[fij]], emodel = emodel)
      sco_X[[fij]] <- PCA_X[[fij]]$PCAscore
    }
    
      Bhat <- est.fun(sco_Y = sco_Y, sco_X = do.call(cbind, sco_X),
                      emodel = emodel, fmodel = fmodel)
      
      Yfit <- matrix(NA, nrow = n, ncol = p)
      for(k in 1:n){
        Xk = do.call(cbind, sco_X)[k,]
        model_k <- pred.fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
        Yfit[k,] <- eval.fd(model_k, seq(gpY[1], gpY[p], length.out = p))
      }
      
      resids <- Y - Yfit
      var_index <- mindex
  }
  
  data <- list()
  data$X <- X
  data$Y <- Y
  
  fpca.results <- list()
  fpca.results$Y <- PCA_Y
  fpca.results$X <- PCA_X
  
  model.details <- list()
  model.details$var.used <- var_index
  model.details$ncompX <- ncompX
  model.details$ncompY <- ncompY
  model.details$nbasisY <- nbasisY
  model.details$nbasisX <- nbasisX
  model.details$gpY <- gpY
  model.details$gpX <- gpX
  model.details$Bhat <- Bhat
  model.details$model <- model
  model.details$emodel <- emodel
  model.details$fmodel <- fmodel
  
  return(list(data = data,
              fitted.values = Yfit,
              residuals = resids,
              fpca.results = fpca.results,
              model.details = model.details))
}
