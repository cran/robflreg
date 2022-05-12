BIC.fun <-
function(Y, Yfit, ncompX, ncompY, emodel){
  n = dim(Y)[1]
  arg_bic = numeric(n)
  for(m in 1:n)
    arg_bic[m] = t(Y[m,] - Yfit[m,]) %*% (Y[m,] - Yfit[m,])
  
  bic_index = sort.int(arg_bic, decreasing = FALSE,
                       index.return = TRUE)
  ntrim = round(0.8 * n)
  index_trunc = bic_index$ix[1:ntrim]
  
  if(emodel == "classical"){
    BIC_val = n * log(sum(arg_bic) / n) +
      (ncompX * ncompY + 1) * log(n)
  }else if(emodel == "robust"){
    BIC_val = ntrim * log(sum(arg_bic[index_trunc]) / ntrim) +
      (ncompX * ncompY + 1) * log(ntrim)
  }
  return(BIC_val)
}
