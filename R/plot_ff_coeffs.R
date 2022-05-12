plot_ff_coeffs <- function(object, b, phi, theta, cex.axis)
{
  coefficients <- object$coefficients
  gpY <- object$gpY
  gpX <- object$gpX
  p <- object$vars
  n <- 1:length(p)
  if(b > length(n))
    stop("Error!! b should be smaller than", " ", length(n), "!")
  
  return(persp3D(z = coefficients[[n[b]]], x = gpX[[n[b]]], y = gpY, ylab="\n Grid point (t)",
                 xlab="\n Grid point (s)", zlab="", main = parse(text = sprintf('beta[%s](s,t)', b)), phi=phi, theta = theta,
                 colkey=FALSE, colvar = coefficients[[n[b]]], ticktype = "detailed",
                 zlim=c(min(coefficients[[n[b]]]),max(coefficients[[n[b]]])), cex.axis=cex.axis))
}
