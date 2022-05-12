generate.sf.data <-
function(n, n.pred, n.gp){
    
    if(!n.pred > 0)
      stop("Error!! Number of predictors must be greater than one!")
    if(!n > 2)
      stop("Error!! Functional variables must have at least two observations!")
    if(!n.gp > 4)
      stop("Error!! The length of grid points must be at least 4 !")
    
    gpX <- seq(0, 1, length.out = n.gp)
    
    cX <- runif(1, min = 1, max = 4)
    
    fX <- fXd <- list()
    for(j in 1:n.pred){
      
      ksi <- list()
      for(i in 1:5){
        ksi[[i]] <- rnorm(n, 1, sd = (cX*i^(-3/2)))
      }
      
      phi <- list()
      for(i in 1:5){
        phi[[i]] <- sin(i * pi * gpX) - cos(i * pi * gpX)
      }
      
      fX[[j]] <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
      fXd[[j]] <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
    }
    
    coef.space <- list()
    coef.space[[1]] <- sin(pi* gpX)
    coef.space[[2]] <- sin(2*pi* gpX)
    coef.space[[3]] <- sin(3*pi* gpX)
    coef.space[[4]] <- sin(4*pi* gpX)
    coef.space[[5]] <- sin(5*pi* gpX)
    coef.space[[6]] <- cos(pi* gpX)
    coef.space[[7]] <- cos(2*pi* gpX)
    coef.space[[8]] <- cos(3*pi* gpX)
    coef.space[[9]] <- cos(4*pi* gpX)
    coef.space[[10]] <- cos(5*pi* gpX)
    
    coef.ind <- numeric()
    if(n.pred <= 10){
      coef.ind <- sample(1:10, n.pred, replace = FALSE)
    }else{
      coef.ind <- sample(1:10, n.pred, replace = TRUE)
    }
    
    vBeta <- coef.space[coef.ind]
    vBeta0 <- vBeta
    
    for(ij in 1:n.pred){
      fX[[ij]] = fdata(fX[[ij]], argvals = gpX)
      vBeta[[ij]] = fdata(runif(1, min = 1, max = 3) * 
                            vBeta[[ij]], argvals = gpX)
    }
    
    err = rnorm(n, mean=0, sd=1)
    
    fY = Reduce("+", lapply(1:length(fX), function(k){inprod.fdata(fX[[k]], vBeta[[k]])}))
    fYe = fY + err
    
    return(list("Y" = fYe, "X" = fXd, "f.coef" = vBeta0))
  }
