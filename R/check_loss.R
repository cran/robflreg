# Check-loss function
check_loss <- function(u, tau)
{
  u * (tau - (u < 0))
}