BIC_fun_JAS <- function(y, yfit, h, nby, nbx, tau)
{
  n = dim(y)[1]
  res = y - yfit
  BIC_val = sum(check_loss(res, tau)) + (nby+nbx+h) * log(n)
  return(BIC_val)
}
