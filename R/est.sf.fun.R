est.sf.fun <-
function(Y, sco_X,
           emodel = c("classical", "robust"),
           fmodel = c("LTS", "MM", "S", "tau")){

    mdl.frame <- data.frame(Y, sco_X)
    for(i in 2:dim(mdl.frame)[2])
      colnames(mdl.frame)[i] = paste("X", i, sep = "")

    emodel <- match.arg(emodel)
    if(emodel != "classical"){
      fmodel <- match.arg(fmodel)
    }else{
      fmodel <- NULL
    }
    if(emodel == "classical"){
      Bhat <- as.matrix(ginv(t(sco_X) %*% sco_X) %*% t(sco_X) %*% Y)
    }else if(emodel == "robust"){
      sco_X <- as.matrix(sco_X)
      Y <- as.matrix(Y)
      if(fmodel == "LTS")
        Bhat <- as.matrix(ltsReg(Y ~ ., data = mdl.frame)$coefficients[-1])
      if(fmodel == "MM")
        Bhat <- as.matrix(lmrob(Y ~ ., data = mdl.frame)$coefficients[-1])
      if(fmodel == "S")
        Bhat <- as.matrix(lmrob.S(x=cbind(1, mdl.frame[,-1]), y = mdl.frame[,1],
                                  control = lmrob.control(nResample = 20), trace.lev=1)$coefficients[-1])
      if(fmodel == "tau")
        Bhat <- as.matrix(FastTau(sco_X, Y)$coefficients)
    }
    return(Bhat)
  }
