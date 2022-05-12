est.fun <-
function(sco_Y, sco_X,
                    emodel = c("classical", "robust"),
                    fmodel = c("MCD", "MLTS", "MM", "S", "tau")){
  
  emodel <- match.arg(emodel)
  if(emodel != "classical"){
    fmodel <- match.arg(fmodel)
  }else{
    fmodel <- NULL
  }
  if(emodel == "classical"){
    Bhat <- ginv(t(sco_X) %*% sco_X) %*% t(sco_X) %*% sco_Y
  }else if(emodel == "robust"){
    sco_X <- as.matrix(sco_X)
    sco_Y <- as.matrix(sco_Y)
    if(fmodel == "MCD")
      fit_model <- mcd(sco_X, sco_Y)
    if(fmodel == "MLTS")
      fit_model <- mlts(sco_X, sco_Y, gamma=.50)
    if(fmodel == "MM")
      fit_model <- MMest_multireg(sco_X, sco_Y, int=F)
    if(fmodel == "S")
      fit_model <- GSest_multireg(sco_X, sco_Y, int=F)
    if(fmodel == "tau")
      fit_model <- taumlmm(X = sco_X, Y = sco_Y, N=dim(sco_X)[1]/2,c1=3, c2=5, ka=1)
    Bhat <- as.matrix(t(fit_model$betaR))
  }
  return(Bhat)
}
