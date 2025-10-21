qcov <- function(y, x, tau, qc.type){
  
  if(qc.type == "dodge"){
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    cf <- numeric()
    for(i in 1:p){
      mod.i <- rq(y~x[,i], tau)
      cf[i] <- mod.i$coefficients[2]
    }
    
    output <- diag(cov(x)) * cf
  }else if(qc.type == "choi"){
    
    y <- as.matrix(y)
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    #qr_slope_y_on_x
    cfx <- numeric()
    for(i in 1:p){
      mod.i <- rq(y~x[,i], tau)
      cfx[i] <- mod.i$coefficients[2]
    }
    
    #qr_slope_x_on_y
    cfy <- numeric()
    for(m in 1:p){
      mod.m <- rq(x[,m]~y, tau)
      cfy[m] <- mod.m$coefficients[2]
    }
    
    # Quantile correlation
    arg1 <- cfy * cfx
    tmp_results <- sqrt(ifelse(arg1>0, arg1, 0))
    tmp_results <- sign(cfy) * tmp_results
    
    # Quantile covariance
    sqrt_var_term <- sqrt(diag(cov(x)) * as.vector(var(y)))
    output <- tmp_results * sqrt_var_term
    
  }else if(qc.type == "li"){
    n <- dim(x)[1]
    psi_value <- tau - (y - (quantile(y, probs=tau) < 0))
    x_center <- scale(x, scale = T)  
    output <- (1.0 / n) * (t(psi_value) %*% x_center)
  }
  return(output)
}