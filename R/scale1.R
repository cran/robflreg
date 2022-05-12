scale1 <-
function(X, k, cgl, tol ,q){
  
  bbb <- rhobi(k, k) / 2
  s0 <- median(abs(X))
  if(q==5){
    s0 <- s0 / 2.1
  }
  if(q==2){
    s0 <- s0 / 1.15 
  }
  if(q==1){
    s0 <- s0 / 0.6745
  }
  s=s0
  if(s0>.00001){   
    gamma <- 1
  }else{
    gamma <- 0
  }
  ii <- 0
  while((gamma > tol)|(ii<100)){
    ii <- ii + 1
    s <- cgl * s0^2 * mean(rhobi(X / s0, k)) / bbb
    s <- sqrt(s)
    gamma <- abs(s - s0) / s0
    s0 <- s
  }
  return(s0)
}
