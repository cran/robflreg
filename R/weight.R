weight <-
function(d, c1, c2, s){
  
  n <- prod(dim(d))
  d <- d / s
  d <- pmax(abs(d), .000000001)
  A2 <- sum(rhobi(d, c2))
  B1 <- sum(psibi(d, c1) * d)
  B2 <- sum(psibi(d, c2) * d)
  CC1 <- (2 * A2 - B2) / n
  CC2 <- B1 / n
  w <- CC1 * psibi(d, c1) + CC2 * psibi(d, c2)
  w <- w / d
  list(w = w)
}
