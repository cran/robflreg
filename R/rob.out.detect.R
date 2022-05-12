rob.out.detect <-
function(object, alpha = 0.01, B = 200, fplot = FALSE){
  
  Y <- object$data$Y
  yarg <- object$model.details$gpY
  residuals <- object$residuals
  dpt <- fdata(residuals, argvals = NULL, rangeval = NULL,
                        names = NULL, fdata2d = FALSE)
  fdepth = depth.mode(dpt)$dep
  ctf <- getCutoff(data = Y, depth = fdepth, alpha = alpha, B = B)
  
  out_indx <- sort(which(fdepth < ctf))
  
  if(fplot){
    yfdt <- fdata(Y, argvals = yarg)
    
    plot(yfdt[-out_indx], lty = 1, ylab = "Response", xlab = "Grid point (t)", 
         main = "", mgp = c(2, 0.5, 0), ylim = range(Y), col = "gray")
    lines(yfdt[out_indx], lty = 1, col = "black")
    legend("topleft", legend = c("Normal", "Outlier"), col = c("grey", "black"), lty = 1, cex = 0.75)
  }
  
  return(cat("outlying functions are:", out_indx))
}
