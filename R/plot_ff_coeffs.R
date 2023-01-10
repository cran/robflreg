plot_ff_coeffs <- function(object, b)
{
  p <- object$vars
  n <- 1:length(p)
  t <- object$gpY
  s <- object$gpX[[n[b]]]
  if(b > length(n))
    stop("Error!! b should be smaller than", " ", length(n), "!")
  coefficients <- object$coefficients[[n[b]]]
  
  
  return(image.plot(t, s, coefficients, main = parse(text = sprintf('hat(beta)[%s](s,t)', b)), 
               xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0)))
}
