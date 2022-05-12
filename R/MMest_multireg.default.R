MMest_multireg.default <-
function(X, Y, int = TRUE, control=MMcontrol(...), na.action=na.omit, ...){
  
  rhobiweight <- function(x,c){
    
    hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
    rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)
    
    return(rho)
  }
  
  scaledpsibiweight <- function(x,c){
    
    hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
    psi <- hulp*(abs(x)<c)
    
    return(psi)
  }
  
  "chi.int" <- function(p, a, c1)
    return(exp(lgamma((p + a)/2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))
  
  
  "loceff.bw" <- function(p, c1){  
    
    alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))   
    beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
    beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
    beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2
    
    return(beta1^2 / alpha1 )
    
  }
  
  csolve.bw.MM <- function(p, eff){
    
    maxit <- 1000
    eps <- 10^(-8)
    diff <- 10^6
    ctest <- -.4024 + 2.2539 * sqrt(p)
    iter <- 1
    while((diff>eps) & (iter<maxit)){
      cold <- ctest
      ctest <- cold * eff / loceff.bw(p,cold)
      diff <- abs(cold-ctest)
      iter <- iter+1
    }
    return(ctest)
    
  }
  
  Y <- as.matrix(Y)
  ynam=colnames(Y)
  q=ncol(na.action(Y))
  if(q < 1L) stop("at least one response needed")
  X <- as.matrix(X)
  xnam=colnames(X)
  if(nrow(Y) != nrow(X))stop("x and y must have the same number of observations")
  YX=na.action(cbind(Y,X))
  Y=YX[,1:q,drop=FALSE]
  X=YX[,-(1:q),drop=FALSE]
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  q <- ncol(Y)
  if(p < 1L) stop("at least one predictor needed")
  if(q < 1L) stop("at least one response needed")
  if(n < (p+q)) stop("For robust multivariate regression the number of observations cannot be smaller than the total number of variables")
  
  interceptdetection <- apply(X==1, 2, all)
  interceptind <- (1:p)[interceptdetection==TRUE]
  if(!any(interceptdetection) & int){
    X <- cbind(rep(1,n),X)
    p <- p + 1    
    interceptind <-1
    if(!is.null(xnam)) colnames(X)[1] <- "(intercept)"
  }
  
  if(is.null(ynam))
    colnames(Y) <- paste("Y",1:q,sep="")
  if(is.null(xnam)){
    colnames(X) <- paste("X",1:p,sep="")
    if(interceptdetection || int){
      colnames(X)[interceptind] <- "(intercept)"
      colnames(X)[-interceptind] <- paste("X",1:(p-1),sep="")
    }  
  }
  eff <- control$eff
  bdp <- control$bdp
  shapeEff <- control$shapeEff
  fastScontrols <- control$fastScontrols
  maxiter <- control$maxIt.MM
  mtol <- control$convTol.MM
  
  c1 <- csolve.bw.MM(m, eff)
  Sresult <- Sest_multireg(X, Y, int=FALSE, bdp, fastScontrols)
  auxscale <- Sresult$scale
  newG <- Sresult$Gamma
  newBeta <- Sresult$coefficients
  newR <- Y - X %*% newBeta
  psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
  newobj <- mean(rhobiweight(psres/auxscale,c1))
  origobj <- newobj
  
  iteration <- 1
  oldobj <- newobj + 1
  while(((oldobj - newobj) > mtol) & (iteration < maxiter)){
    oldobj <- newobj
    w <- scaledpsibiweight(psres/auxscale,c1)
    wbig <- matrix(rep(w,p),ncol=p)	
    wX <- X * wbig	
    newBeta <- ginv(crossprod(wX, X)) %*% crossprod(wX, Y)
    newG <- cov.wt(newR, wt=w, center=FALSE)$cov
    newG <- det(newG)^(-1/m)*newG
    newR <- Y - X %*% newBeta
    psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
    newobj <- mean(rhobiweight(psres/auxscale,c1))
    iteration <- iteration+1
  }
  
  if(newobj <= origobj){
    resBeta <- newBeta
    resshape <- newG
    rescovariance <- newG*auxscale^2
  }else{
    resBeta <- Sresult$coefficients
    resshape <- Sresult$Gamma
    rescovariance <- Sresult$Gamma*auxscale^2
  }
  
  return(list(betaR=t(resBeta)))       
}
