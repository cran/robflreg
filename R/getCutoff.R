getCutoff <-
function(data, depth, alpha, B){
  
  Call = numeric()
  for(nb in 1:B){
    prbs = depth/sum(depth)
    indx = sample(1:length(depth), length(depth), replace = TRUE, prob = prbs)
    
    Yboot = data[indx,] + 0.05 * matrix(rmvnorm(dim(data)[1], mean = rep(0, dim(data)[2]), 
                                                sigma = cov(data)), ncol = dim(data)[2])
    boot_depths = depth.mode(fdata(Yboot, argvals = NULL, rangeval = NULL,
                                   names = NULL, fdata2d = FALSE))$dep
    Call[nb] = quantile(boot_depths, probs = alpha)
  }
  
  return(median(Call))
}
