rob.sf.reg <-
  function(Y, X, X.scl = NULL, emodel = c("classical", "robust"), fmodel = c("LTS", "MM", "S", "tau"),
           nbasis = NULL, gp = NULL, ncomp = NULL){

    if(!is.null(X.scl)){
      X.scl <- X.scl0 <- as.matrix(X.scl)
    }else{
      X.scl0 <- NULL
    }
    if(!is.matrix(Y))
      Y <- as.matrix(Y)
    if(!is.list(X))
      stop("Error!! X must be a list!")
    if(sum(sapply(1:length(X),function(i){!is.matrix(X[[i]])})))
      stop("Error!! X must be a list and all of its components must be matrix!")
    dimX <- sapply(1:length(X),function(i){dim(X[[i]])[1]})
    if((length(unique(dimX))!=1))
      stop("Error!! All the components of X must be matrix and they must have the same number of rows!")
    if(dimX[1] != dim(Y)[1])
      stop("Error!! Response and functional predictor variables must have the same number of rows!")
    if(!is.null(X.scl)){
      if(dim(X.scl)[1] != dim(Y)[1])
        stop("Error!! Response and scalar predictor variables must have the same number of rows!")
    }
    if(!is.null(gp)){
      if(!is.list(gp))
        stop("Error!! gp must be a list!")
      lengpX <- sapply(1:length(gp), function(i){length(gp[[i]])})
      nc <- sapply(1:length(X), function(i){dim(X[[i]])[2]})
      if(length(gp) != length(X))
        stop("Error!! The lengths of X and gp must be equal!")
      if(!all(lengpX == nc))
        stop("Error!! The number of columns of each component of X must be equal to the length of corresponding component of gp!")
    }
    if(!is.null(nbasis)){
      if(!is.numeric(nbasis))
        stop("Error!! nbasis must be a numeric vector!")
      if(!all(nbasis >= 4))
        stop("Error!! Each component of nbasis must be greater than three to apply cubic B-spline basis expansion!")
      if(length(nbasis) != length(X))
        stop("Error!! The lengths of X and nbasis must be equal!")
    }
    if(!is.null(ncomp)){
      if(!is.numeric(ncomp))
        stop("Error!! ncomp must be a numeric vector!")
      if(!all(ncomp > 0))
        stop("Error!! Each component of ncomp (the numbers of principal components) must be greater than or equal to one!")
      if(length(ncomp) != length(X))
        stop("Error!! The lengths of X and ncomp must be equal!")
    }
    if(!emodel %in% c("classical", "robust"))
      stop("Error!! emodel must be one of the followings: classical or robust!")
    if(emodel == "robust"){
      if(!fmodel %in% c("LTS", "MM", "S", "tau"))
        stop("Error!! model must be one of the followings: LTS, MM, S, and tau!")
    }

    emodel <- match.arg(emodel)
    if(emodel != "classical"){
      fmodel <- match.arg(fmodel)
    }else{
      fmodel <- NULL
    }

    np <- length(X)
    n <- dim(Y)[1]
    px <- numeric()
    for(ip in 1:np)
      px[ip] <- dim(X[[ip]])[2]

    if(is.null(nbasis)) nbasis = rep(min(20, n/4), np)
    if(is.null(gp)){
      gp <- list()
      for(iip in 1:np)
        gp[[iip]] = seq(0, 1, length.out = px[iip])
    }

    if(emodel == "classical"){
      mean_Y <- mean(Y)
      Y <- scale(Y, scale = F)
    }
    if(emodel == "robust"){
      mean_Y <- median(Y)
      Y <- scale(Y, center = mean_Y, scale = F)
    }
    if(!is.null(X.scl)){
      if(emodel == "classical")
        mean_X.scl <- apply(X.scl, 2, mean)
      if(emodel == "robust")
        mean_X.scl <- apply(X.scl, 2, median)
      X.scl <- scale(X.scl, center = mean_X.scl, scale = F)
    }

    PCA_X <- list()
    sco_X <- list()
    ncompXv <- numeric()
    for(fij in 1:np){
      PCA_X[[fij]] <- getPCA(data = X[[fij]], nbasis = nbasis[fij], ncomp = ncomp[fij],
                             gp = gp[[fij]], emodel = emodel)
      sco_X[[fij]] <- PCA_X[[fij]]$PCAscore
      ncompXv[fij] <- PCA_X[[fij]]$ncomp
    }

    if(is.null(ncomp))
      ncomp <- ncompXv

    X.preds <- cbind(X.scl, do.call(cbind, sco_X))

    Bhat <- est.sf.fun(Y = Y, sco_X = X.preds,
                       emodel = emodel, fmodel = fmodel)

    Yfit <- X.preds %*% Bhat + mean_Y
    resids <- Y - Yfit

    if(!is.null(X.scl)){
      Bshat <- Bhat[1:dim(X.scl)[2],]
      Bhat <- Bhat[-(1:dim(X.scl)[2])]
    }else{
      Bshat <- NULL
      Bhat <- Bhat
    }

    data <- list()
    data$X <- X
    data$Y <- Y
    data$X.scl <- X.scl0

    fpca.results <- PCA_X

    model.details <- list()
    model.details$mean_Y <- mean_Y
    model.details$ncomp <- ncomp
    model.details$nbasis <- nbasis
    model.details$gp <- gp
    model.details$Bhat <- Bhat
    model.details$Bshat <- Bshat
    model.details$emodel <- emodel
    model.details$fmodel <- fmodel

    return(list(data = data,
                fitted.values = Yfit,
                residuals = resids,
                fpca.results = fpca.results,
                model.details = model.details))
  }
