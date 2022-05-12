tau <-
function(X, c1, c2, ka, cgl, tol, scal){
  
  if(scal==-1){
    s <- scale1(X, c1, cgl, tol, 1)
  }else{
    s <- scal
  }
  t <- ( 1 / (6 * ka)) * (s^2) * (c2^2) * mean(rhobi(X / s, c2))
  t <- sqrt(t)
  list(t=t, s=s)
}
