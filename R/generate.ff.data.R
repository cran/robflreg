generate.ff.data <-
function(n.pred, n.curve, n.gp, out.p = 0){
  
  if(!n.pred > 0)
    stop("Error!! Number of predictors must be greater than one!")
  if(!n.curve > 2)
    stop("Error!! Functional variables must have at least two observations!")
  if(!n.gp > 4)
    stop("Error!! The length of grid points must be at least 4 !")
  if(out.p < 0 | out.p > 1)
    stop("Error!! Outlier percentage must be between 0 and 1!")
  
  gpX <- gpY <- seq(0, 1, length.out = n.gp)
  
  cX <- runif(1, min = 1, max = 4)
  
  fX <- list()
  for(j in 1:n.pred){
    
    ksi <- list()
    for(i in 1:5){
      ksi[[i]] <- rnorm(n.curve, 1, sd = (cX*i^(-1/2)))
    }
    
    phi <- list()
    for(i in 1:5){
      phi[[i]] <- sin(i * pi * gpX) - cos(i * pi * gpX)
    }
    
    fX[[j]] <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  }
  
  coef.space <- list()
  coef.space[[1]] <- function(s,t) sin(2 * pi * s) * sin(pi * t)
  coef.space[[2]] <- function(s,t) cos(3 * pi * s) * sin(2 * pi * t)
  coef.space[[3]] <- function(s,t) sin(pi * s) * cos(0.5 * pi * t)
  coef.space[[4]] <- function(s,t) cos(3/2 * pi * s) * cos(3/2 * pi * t)
  coef.space[[5]] <- function(s,t) cos(pi * s) * cos(6 * pi * t)
  coef.space[[6]] <- function(s,t) exp(-(s - 0.5)^2) * exp(-2 * (t - 1)^2)
  coef.space[[7]] <- function(s,t) exp(-3*(s - 0.5)^2) * exp(-4 * (t - 1)^2)
  coef.space[[8]] <- function(s,t) (s - 0.5)^2 * (t - 0.5)^2
  coef.space[[9]] <- function(s,t) (s - 0.5)^2 * sin(2 * pi * t)
  coef.space[[10]] <- function(s,t) sqrt(s) * sqrt(2 * t)
  
  fBeta <- list()
  coef.ind <- numeric()
  if(n.pred <= 10){
    coef.ind <- sample(1:10, n.pred, replace = FALSE)
  }else{
    coef.ind <- sample(1:10, n.pred, replace = TRUE)
  }
  
  fBeta <- coef.space[coef.ind]
  vBeta <- lapply(1:n.pred, function(k){outer(gpX, gpY, fBeta[[k]])})
  for(iv in 1:n.pred)
    vBeta[[iv]] <- vBeta[[iv]] * runif(1, min = 1, max = 3)

  fY <- Reduce("+", lapply(1:n.pred, function(k){fX[[k]] %*% vBeta[[k]] / n.gp})) +
    r_ou(n=n.curve, t = gpY, mu=0, alpha = 1, sigma = 1,
         x0=rnorm(n=n.curve, mean = 0, sd = 1/sqrt(2)))$data
  
  out.indx <- NULL
  
  if(out.p > 0){

    fX.out <- list()
    for(j in 1:n.pred){
      
      ksi.out <- list()
      for(i in 1:5){
        ksi.out[[i]] <- rnorm(n.curve, 1, sd = (cX*i^(-3/2)))
      }
      
      phi.out <- list()
      for(i in 1:5){
        phi.out[[i]] <- 2*sin(i * pi * gpX) - cos(i * pi * gpX)
      }
      
      fX.out[[j]] <- Reduce("+", lapply(1:5, function(k){ksi.out[[k]] %*% t(phi.out[[k]])}))
    }
    
    fBeta.out <- list()
    coef.ind.out <- numeric()
    if(n.pred <= 10){
      coef.ind.out <- sample(1:10, n.pred, replace = FALSE)
    }else{
      coef.ind.out <- sample(1:10, n.pred, replace = TRUE)
    }
    
    fBeta.out <- coef.space[coef.ind.out]
    vBeta.out <- lapply(1:n.pred, function(k){outer(gpX, gpY, fBeta.out[[k]])})
    for(iv in 1:n.pred)
      vBeta.out[[iv]] <- vBeta.out[[iv]] * runif(1, min = 1, max = 2)
    
    fY.out <- Reduce("+", lapply(1:n.pred, function(k){fX.out[[k]] %*% vBeta.out[[k]] / n.gp})) +
      r_ou(n=n.curve, t = gpY, mu=0, alpha = 1, sigma = 1,
           x0=rnorm(n=n.curve, mean = 0, sd = 1/sqrt(2)))$data
    
    
    nout <- round(n.curve * out.p)
    out.indx <- sample(1:n.curve, nout)
    
    fY[out.indx,] <- fY.out[out.indx,]
    for(io in 1:n.pred)
      fX[[io]][out.indx,] <- fX.out[[io]][out.indx,]
  }

  return(list("Y" = fY, "X" = fX, "f.coef" = vBeta, out.indx = out.indx))
  
}
