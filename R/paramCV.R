paramCV <- function(y, x, hs, tau, nbys, nbxs, gpy, gpx,
                    nfold, qc.type)
{

  n <- dim(x)[1]
  p <- dim(x)[2]

  folds <- cvFolds(n, K = nfold, R = 1, type = "random")
  CVmat <- as.matrix(expand.grid(hs, nbys, nbxs))
  CVmat <- cbind(CVmat, NA)
  colnames(CVmat) <- c("hs", "nbys", "nbxs", "BIC")
  CVmat <- CVmat[CVmat[,2] >= CVmat[,1],]

  data <- cbind(y, x)

  for(i in 1:dim(CVmat)[1]){
    h <- CVmat[,1][i]
    nby <- CVmat[,2][i]
    nbx <- CVmat[,3][i]
    BIC_mat <- numeric()

    for(f in 1:nfold){
      try({
        dtrain <- data[folds$which!=f,]
        dtest <- data[folds$which==f,]

        yf <- dtrain[,(1:dim(y)[2])]
        xf <- dtrain[,-(1:dim(y)[2])]
        yft <- dtest[,(1:dim(y)[2])]
        xft <- dtest[,-(1:dim(y)[2])]

        BS.sol.x <- getAmat(data = xf, nbf = nbx, gp = gpx)
        x1 <- BS.sol.x$Amat
        BS.sol.y <- getAmat(data = yf, nbf = nby, gp = gpy)
        y1 <- BS.sol.y$Amat

        trainmod <- pqr(y = y1, x = x1, tau = tau, h = h, qc.type = qc.type)
        Bs.sol.xnew <- getAmat(data = xft, nbf = nbx, gp = gpx)
        x1new <- Bs.sol.xnew$Amat

        predictions <- (cbind(1,x1new) %*% trainmod$d.coef) %*% solve(BS.sol.y$sinp_mat) %*% t(BS.sol.y$evalbase)
        BIC_mat[f] <- BIC_fun_JAS(smooth_fun(yft), predictions, h, nby, nbx, tau)
      }, silent = TRUE)
    }
    CVmat[,4][i] <- mean(BIC_mat)
  }

  optima <- CVmat[which.min(CVmat[,4]),][1:3]
  optimh <- optima[1]
  optimnby <- optima[2]
  optimnbx <- optima[3]
  return(list(h = optimh, nby = optimnby, nbx = optimnbx))
}
