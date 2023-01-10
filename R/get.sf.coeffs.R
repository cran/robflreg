get.sf.coeffs <-
function(object){

    np <- length(object$data$X)
    ncX <- object$model.details$ncomp
    gp <- object$model.details$gp
    Bhat <- object$model.details$Bhat
    Bshat <- object$model.details$Bshat
    evalx <- list()
    fcomp_X <- list()
    for(i in 1:np){
      evalx[[i]] <- object$fpca.results[[i]]$evalbase
      fcomp_X[[i]] <- object$fpca.results[[i]]$PCAcoef$coefs
    }

    coeffs = list()
    km = 1
    for(j in 1:np){
      coeffs[[j]] = evalx[[j]] %*% (fcomp_X[[j]] %*% Bhat[km: (km+ncX[j]-1),])
      km = cumsum(ncX)[j]+1
    }

    return(list(gp = gp, coefficients = coeffs,
                scl.coefficients = Bshat))
  }
