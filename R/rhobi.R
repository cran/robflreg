rhobi <-
function(u, ccc){
  
  w <- abs(u) <= ccc
  v <- (u^2 / (2) * (1 - (u^2 / (ccc^2)) + (u^4 / (3 * ccc^4)))) * w + (1 - w) * (ccc^2 / 6)
  v <- v * 6 / ccc^2
  return(v)
}
