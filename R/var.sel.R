var.sel <-
function(Y, X, nbasisY, nbasisX, ncompX = NULL, ncompY = NULL, gpY, gpX, 
                    emodel = c("classical", "robust"), fmodel = c("MCD", "MLTS", "MM", "S", "tau")){
  emodel <- match.arg(emodel)
  if(emodel != "classical"){
    fmodel <- match.arg(fmodel)
  }else{
    fmodel <- NULL
  }
  np <- length(X)
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  PCA_Y <- getPCA(data = Y, nbasis = nbasisY, ncomp = ncompY,
                       gp = gpY, emodel = emodel)
  
  sco_Y <- PCA_Y$PCAscore
  comp_Y <- PCA_Y$PCAcoef
  mean_Y <- PCA_Y$meanScore
  ncompYv <- PCA_Y$ncomp
  
  sco_X <- list()
  ncompXv <- numeric()
  for(ij in 1:np){
    PCA_X <- getPCA(data = X[[ij]], nbasis = nbasisX[ij], ncomp = ncompX,
                        gp = gpX[[ij]], emodel = emodel)
    sco_X[[ij]] <- PCA_X$PCAscore
    ncompXv[ij] <- PCA_X$ncomp
  }

  BIC_individuals <- numeric()
  
  for(ind in 1:np){
    Bhat <- est.fun(sco_Y = sco_Y, sco_X = sco_X[[ind]],
                   emodel = emodel, fmodel = fmodel)
    
    Yhat <- matrix(, nrow = n, ncol = p)
    
    for(k in 1:n){
      Xk = sco_X[[ind]][k,]
      model_k <- pred.fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
      Yhat[k,] <- eval.fd(model_k, seq(0, 1, length.out = p))
    }
    
    BIC_individuals[ind] <- BIC.fun(Y = Y, Yfit = Yhat, ncompX = ncompXv,
                                   ncompY = ncompYv, emodel = emodel)
  }
  BIC_order <- order(BIC_individuals)
  
  main_model_start <- sco_X[[BIC_order[1]]]
  BIC_forw <- min(BIC_individuals)
  
  X_next <- c(BIC_order[1], rep(NA, (length(BIC_order)-1)))
  X_out <- c(which.min(BIC_individuals))
  
  for(f1 in 2:np){
    BIC_sel <- rbind(subset(BIC_order, !(BIC_order %in% X_out)), NA)
    for(f2 in 1:ncol(BIC_sel)){
      sco_X_forw <- cbind(main_model_start, sco_X[[BIC_sel[1,f2]]])
      
      Bhat <- est.fun(sco_Y = sco_Y, sco_X = sco_X_forw,
                     emodel = emodel, fmodel = fmodel)
      
      Yhat <- matrix(, nrow = n, ncol = p)
      
      for(k in 1:n){
        Xk <- sco_X_forw[k,]
        model_k <- pred.fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
        Yhat[k,] <- eval.fd(model_k, seq(0, 1, length.out = p))
      }
      
      BIC_sel[2,f2] <- BIC.fun(Y = Y, Yfit = Yhat, ncompX = ncompXv,
                              ncompY = ncompYv, emodel = emodel)
    }
    
    BIC_next <- BIC_sel[2,][which.min(BIC_sel[2,])]
    BIC_next2 <- BIC_sel[1,][which.min(BIC_sel[2,])]
    
    X_out <- c(X_out, BIC_next2)
    
    if(BIC_next/BIC_forw < 0.98){
      main_model_start <- cbind(main_model_start, sco_X[[BIC_next2]])
      BIC_forw <- BIC_next
      X_next[f1] <- BIC_next2
    }else{
      main_model_start <- main_model_start
      BIC_forw <- BIC_forw
      X_next[f1] <- X_next[f1]
    }
  }
  
  selected_main <- sort(subset(X_next, !(X_next %in% NA)))

  return(list(maine = selected_main))
}
