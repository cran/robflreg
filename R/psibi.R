psibi <-
function(u, ccc){
  
  w <- abs(u) <= ccc
  v <- w * u * ((1 - (u / ccc)^2)^2)
  v <- v * 6 / ccc^2
  return(v)
}
