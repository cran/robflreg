# Function to smooth raw data
smooth_fun <- function(data, argvals = NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if(is.null(argvals))
    argvals <- seq(0, 1, length.out = p)
  nbasis <- min(40, round(p/5))
  basis.obs <- create.bspline.basis(range(argvals), nbasis)
  
  sdata <- matrix(, nrow = n, ncol = p)
  for(i in 1:n){
    xs <- smooth.basis(argvals = argvals, y= c(data[i,]), fdParobj = basis.obs)
    xfd <- xs$fd
    sdata[i,] <- eval.fd(argvals, xfd)
  }
  
  return(sdata)
}
