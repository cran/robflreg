plot_sf_coeffs <-
function(object, b)
{
  coefficients <- object$coefficients
  gp <- object$gp
  n <- 1:length(coefficients)
  if(b > length(n))
    stop("Error!! b should be smaller than", " ", length(n), "!")
  
  return(plot(gp[[b]], coefficients[[b]], type = "l", xlab = "\n Grid point (s)",
       ylab = "", main = bquote(hat(beta)[.(b)](s))))
}

