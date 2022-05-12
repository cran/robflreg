pred.fun <-
function(comp_Y, sco_X, Bhat){
  
  ncomp <- dim(comp_Y$coefs)[2]
  nest <- t(sco_X %*% Bhat)
  if(ncomp == 1){
    nh <- nest[1] * comp_Y[1,]
  }else{
    nh <- nest[1] * comp_Y[1,]
    for(j in 2:ncomp){
      nh <- nh + nest[j] * comp_Y[j,]  
    }
  }
  return(nh)
}
