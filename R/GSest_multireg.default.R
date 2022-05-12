GSest_multireg.default <-
function(X,Y, int = TRUE, bdp=.5, control=GScontrol(...), na.action=na.omit, ...){

  resmultiGS<-function(x,y,initialbeta,initialgamma,k,b,cc,initialscale,convTol){  

    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)
    
    beta <- initialbeta
    res <- y-x%*%beta
    places <- t(combn(1:n,2))
    ngroot <- nrow(places)
    term1resvector <- as.matrix(res[places[,1],])
    term2resvector <- as.matrix(res[places[,2],])
    diffres <- term1resvector-term2resvector
    rdis <- sqrt(mahalanobis(diffres, rep(0,q), initialgamma))
    term1xvector <- x[places[,1],]
    term2xvector <- x[places[,2],]
    diffx <- term1xvector-term2xvector
    term1yvector <- y[places[,1],]
    term2yvector <- y[places[,2],]
    diffy <- term1yvector-term2yvector
    if(initialscale > 0)
      scale <- initialscale
    else
      scale <- median(rdis)/.6745
    
    gamma <- initialgamma
    
    for(i in 1:k){
      scale <- sqrt( scale^2 * mean( rhobiweight( rdis / scale, cc ) ) / b)
      
      w <- scaledpsibiweight( rdis/scale, cc )
      
      if(p>0){
        wbig <- matrix(rep(w,p),ncol=p)	
        wdiffx <- diffx * wbig
        newbeta <- ginv(crossprod(wdiffx, diffx)) %*% crossprod(wdiffx, diffy)
      }else
        newbeta <- beta
      
      newgamma <- cov.wt(diffres, wt=w, center=FALSE)$cov
      newgamma <- det(newgamma)^(-1/q)*newgamma
      
      if(p>0){  
        if(sum((newbeta-beta)^2)/sum(beta^2) < convTol)       break
      }else{
        if(sum(sum((newgamma-gamma)^2))/sum(sum(gamma^2)) < convTol)       break
      }
      res <- y-x%*%newbeta 
      term1resvector <- as.matrix(res[places[,1],])
      term2resvector <- as.matrix(res[places[,2],])
      diffres <- term1resvector-term2resvector
      rdis <- sqrt(mahalanobis(diffres, rep(0,q), newgamma))
      gamma <- newgamma
      beta <- newbeta
    }    
    
    return(list( betarw = beta, gammarw = newgamma, scalerw = scale,rdis=rdis ))
  }
  
  scalemultiGS <- function(u, b, cc, initialsc){

    if(initialsc==0){
      initialsc <- median(abs(u))/.6745}
    
    maxit <- 200
    sc <- initialsc
    i <- 0
    eps <- 1e-20
    err <- 1
    while((i < maxit ) & (err > eps)){
      sc2 <- sqrt( sc^2 * mean( rhobiweight( u / sc, cc ) ) / b)
      err <- abs(sc2/sc - 1)
      sc <- sc2
      i <- i+1
    }
    
    return(sc)
  }

  rhobiweight <- function(x,c){

    hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
    rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)
    
    return(rho)
  }
  

  psibiweight <- function(x,c){

    hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
    psi <- hulp*(abs(x)<c)
    
    return(psi)
  }
  

  scaledpsibiweight <- function(x,c){

    hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
    psi <- hulp*(abs(x)<c)
    
    return(psi)
  }
  
  vecop <- function(mat){

    nr <- nrow(mat)
    nc <- ncol(mat)
    
    vecmat <- rep(0,nr*nc)
    for(col in 1:nc){
      startindex <- (col-1)*nr+1
      vecmat[startindex:(startindex+nr-1)] <- mat[,col]
    }
    return(vecmat)
  }
  

  reconvec <- function(vec,ncol){

    lcol <- length(vec)/ncol
    rec <- matrix(0,lcol,ncol)
    for(i in 1:ncol)
      rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]
    
    return(rec)
  }
  

  erf <- function(x){
    uitk <- 2*pnorm(x*sqrt(2))-1 
    return(uitk)
  }
  
  TbscGS <- function(alpha,p){

    talpha <- sqrt(qchisq(1-alpha,p))
    maxit <- 1000 
    eps <- 1e-8
    diff <- 1e6
    ctest <- talpha
    iter <- 1
    while((diff>eps) * (iter<maxit)){
      cold <- ctest
    if(alpha >= 0.50){
      ctest <- TbsbGS(cold,p)/(1-alpha^2)}
    else{
      ctest <- TbsbGS(cold,p)/(1-(1-alpha)^2)}
    diff <- abs(cold-ctest)
    iter <- iter+1
    }
    return(ctest)
  }
  
  TbsbGS <- function(c,p){
    
    if(p==1){ 
      y1 <- -1/3*(60*exp(-1/4*c^2)*c+exp(-1/4*c^2)*c^5-60*pi^(1/2)*erf(1/2*c)-3*pi^(1/2)*
                     erf(1/2*c)*c^4-8*exp(-1/4*c^2)*c^3+18*c^2*pi^(1/2)*erf(1/2*c))/(c^4*pi^(1/2))
    y2 <- -1/6*c^2*erf(1/2*c)+1/6*c^2
    res <- (6/c)*(y1+y2)
    }else{
      if(p==2){
      tus <- -1/6*(-384+96*c^2-12*c^4+384*exp(-1/4*c^2)+exp(-1/4*c^2)*c^6)/c^4+1/6*c^2*exp(-1/4*c^2)
      res <- (6/c)*tus
    }else{
      if(p==3){
        tus <- -1/6*(2*exp(-1/4*c^2)*c^5*pi-40*exp(-1/4*c^2)*c^3*pi+840*exp(-1/4*c^2)*c*pi-840*erf(1/2*c)*pi^(3/2)+180*erf(1/2*c)*pi^(3/2)*c^2+exp(-1/4*c^2)*c^7*pi-18*erf(1/2*c)*pi^(3/2)*c^4)/(pi^(3/2)*c^4)+1/6*c^2*(exp(-1/4*c^2)*c*pi^2-pi^(5/2)*erf(1/2*c)+pi^(5/2))/pi^(5/2)
      res <- (6/c)*tus
      }else{
        if(p==4){
          tus <- -1/24*(-96*c^4-6144+1152*c^2+6144*exp(-1/4*c^2)+384*c^2*exp(-1/4*c^2)+4*exp(-1/4*c^2)*c^6+exp(-1/4*c^2)*c^8)/c^4+1/24*c^2*exp(-1/4*c^2)*(4+c^2);
        res <- (6/c)*tus
        }else{
          if(p==5){
            tus <- -1/36*(12*exp(-1/4*c^2)*c^5*pi^2+exp(-1/4*c^2)*c^9*pi^2+2520*erf(1/2*c)*pi^(5/2)*c^2+6*exp(-1/4*c^2)*c^7*pi^2-180*erf(1/2*c)*pi^(5/2)*c^4-15120*pi^(5/2)*erf(1/2*c)+15120*exp(-1/4*c^2)*c*pi^2)/(pi^(5/2)*c^4)+1/36*c^2*(pi^4*exp(-1/4*c^2)*c^3+6*pi^4*exp(-1/4*c^2)*c-6*pi^(9/2)*erf(1/2*c)+6*pi^(9/2))/pi^(9/2)
          res <- (6/c)*tus
          }else{
            if(p==6){
              tus <- -1/192*(-1152*c^4-122880+18432*c^2+384*exp(-1/4*c^2)*c^4+122880*exp(-1/4*c^2)+12288*exp(-1/4*c^2)*c^2+32*exp(-1/4*c^2)*c^6+8*exp(-1/4*c^2)*c^8+exp(-1/4*c^2)*c^10)/c^4+1/192*c^2*exp(-1/4*c^2)*(32+8*c^2+c^4)
            res <- (6/c)*tus
            }else{
              if(p==7){
                tus <- -1/360*(60*exp(-1/4*c^2)*c^7*pi^3+45360*erf(1/2*c)*pi^(7/2)*c^2+504*exp(-1/4*c^2)*c^5*pi^3+10*exp(-1/4*c^2)*c^9*pi^3+332640*exp(-1/4*c^2)*c*pi^3+10080*exp(-1/4*c^2)*c^3*pi^3-2520*erf(1/2*c)*pi^(7/2)*c^4+exp(-1/4*c^2)*c^11*pi^3-332640*erf(1/2*c)*pi^(7/2))/(pi^(7/2)*c^4)+1/360*c^2*(pi^6*exp(-1/4*c^2)*c^5+10*pi^6*exp(-1/4*c^2)*c^3+60*pi^6*exp(-1/4*c^2)*c-60*pi^(13/2)*erf(1/2*c)+60*pi^(13/2))/pi^(13/2)
              res <- (6/c)*tus
              }else{
                if(p==8){
                  tus <- -1/2304*(-18432*c^4-2949120+368640*c^2+18432*exp(-1/4*c^2)*c^4+2949120*exp(-1/4*c^2)+368640*exp(-1/4*c^2)*c^2+768*exp(-1/4*c^2)*c^6+96*exp(-1/4*c^2)*c^8+exp(-1/4*c^2)*c^12+12*exp(-1/4*c^2)*c^10)/c^4+1/2304*c^2*exp(-1/4*c^2)*(384+96*c^2+12*c^4+c^6)
                res <- (6/c)*tus
                }else{
                  if(p==9){
                    tus <- -1/5040*(23184*exp(-1/4*c^2)*c^5*pi^4+exp(-1/4*c^2)*c^13*pi^4+140*exp(-1/4*c^2)*c^9*pi^4+14*exp(-1/4*c^2)*c^11*pi^4+997920*erf(1/2*c)*pi^(9/2)*c^2-8648640*erf(1/2*c)*pi^(9/2)+1224*exp(-1/4*c^2)*c^7*pi^4+443520*exp(-1/4*c^2)*c^3*pi^4-45360*erf(1/2*c)*pi^(9/2)*c^4+8648640*exp(-1/4*c^2)*c*pi^4)/(pi^(9/2)*c^4)+1/5040*c^2*(pi^8*exp(-1/4*c^2)*c^7+14*pi^8*exp(-1/4*c^2)*c^5+140*pi^8*exp(-1/4*c^2)*c^3+840*pi^8*exp(-1/4*c^2)*c-840*pi^(17/2)*erf(1/2*c)+840*pi^(17/2))/pi^(17/2)
                  res <-(6/c)*tus
                  }else{
                    if(p==10){
                      tus <- -1/36864*(-368640*c^4+8847360*c^2-82575360+737280*c^4*exp(-1/4*c^2)+11796480*c^2*exp(-1/4*c^2)+82575360*exp(-1/4*c^2)+30720*exp(-1/4*c^2)*c^6+1920*exp(-1/4*c^2)*c^8+16*exp(-1/4*c^2)*c^12+192*exp(-1/4*c^2)*c^10+exp(-1/4*c^2)*c^14)/c^4+1/36864*c^2*exp(-1/4*c^2)*(6144+1536*c^2+192*c^4+16*c^6+c^8)
                    res <- (6/c)*tus
                    }else{
                      if(p>10){
                        res <- TbsbGSbenader(c,p)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
  
  
  TbsbGSbenader <- function(c,p){
    
    integralpart1 <- function(r){((r^2)/2-(r^4)/(2*c^2)+(r^6)/(6*c^4))*exp((-r^2)/4)*r^(p-1)}
    
    integralpart2<- function(r){exp((-r^2)/4)*r^(p-1)}
    
    tus1<-(2/(gamma(p/2)*2^p))*integrate(integralpart1,0,c)$value
    tus2<-((c^2)/(gamma(p/2)*3*2^p))*integrate(integralpart2,c,Inf)$value
    tustot<-tus1+tus2
    res<-(6/c)*tustot
    return(res)
  }
  
  IRLSlocation <- function(xmat,covmat,bdp,cc){
    
    xmat <- as.matrix(xmat)
    n <- nrow(xmat)
    p <- ncol(xmat)
    neem <- sample(n,p+1)
    xsub <- as.matrix(xmat[neem,])
    
    initmu <- apply(xsub,2,mean)
    initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))
    
    initobj <- mean(rhobiweight(initrdis,cc))
    
    weights <- scaledpsibiweight(initrdis,cc)
    
    itertest <- 0
    while((sum(weights)==0) && (itertest<500)){
      
      neem <- sample(n,p+1)
      xsub <- as.matrix(xmat[neem,])
      
      initmu <- apply(xsub,2,mean)
      initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))
      initobj <- mean(rhobiweight(initrdis,cc))
      weights <- scaledpsibiweight(initrdis,cc)
      itertest <- itertest + 1
    }
    
    if(itertest==500) stop("could not find suitable starting point for IRLS for intercept")   
    
    sqrtweights <- sqrt(weights)
    munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))
    
    rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))
    
    objnieuw <- mean(rhobiweight(rdisnieuw,cc))
    
    iter <- 0
    while(((abs(initobj/objnieuw)-1) > 10^(-15)) && (iter < 100)){
      initobj <- objnieuw
      initmu <- munieuw
      initrdis <- rdisnieuw
      weights <- scaledpsibiweight(initrdis,cc)
      sqrtweights <- sqrt(weights)
      munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))
      rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))
      objnieuw <- mean(rhobiweight(rdisnieuw,cc))
      iter <- iter + 1
    }
    return(munieuw)
  }

  Y <- as.matrix(Y)
  ynam=colnames(Y)
  q=ncol(na.action(Y))
  if(q < 1L) stop("at least one response needed")
  X <- as.matrix(X)
  xnam=colnames(X)
  if(nrow(Y) != nrow(X))stop("x and y must have the same number of observations")
  YX = na.action(cbind(Y,X))
  Y=YX[,1:q,drop=FALSE]
  X=YX[,-(1:q),drop=FALSE]
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  if((p < 1L) && !int) stop("at least one predictor needed")
  if(q < 1L) stop("at least one response needed")
  if(n < (p+q)) stop("For robust multivariate regression the number of observations cannot be smaller than the total number of variables")
  
  interceptdetection <- apply(X==1, 2, all)
  if(any(interceptdetection)) int=TRUE
  zonderint <- (1:p)[interceptdetection==FALSE]
  xzonderint <- X[,zonderint,drop=FALSE]
  X <- xzonderint
  p<-ncol(X)

  if(is.null(ynam))
    colnames(Y) <- paste("Y",1:q,sep="")
  if(is.null(xnam)&& (p>=1L)) 
    colnames(X) <- paste("X",1:p,sep="")
  
  ngroot <- choose(n,2)
  cc <- TbscGS(bdp,q)
  b <- (cc/6)*TbsbGS(cc,q)
  nsamp <- control$nsamp
  bestr <- control$bestr
  k <- control$k
  convTol <- control$convTol
  maxIt <- control$maxIt
  
  bestbetas <- matrix(0,p*q,bestr)
  bestgammas <- matrix(0,q*q,bestr)
  bestscales <- 1e20 * rep(1,bestr)
  sworst <- 1e20
  xextra <- cbind(rep(1,n),X)
  for(i in 1:nsamp){ 
    
    rankR <- 0
    itertest <- 0
    
    while((rankR < q) && (itertest<200)){
      ranset <- sample(n,p+q+1)
      xj <- xextra[ranset,,drop=FALSE]
      yj <- Y[ranset,,drop=FALSE]
      beta <- ginv(crossprod(xj)) %*% crossprod(xj,yj)
      res <- yj - xj %*% beta
      qrRj <- qr(res)
      rankR <- qrRj$rank
      itertest <- itertest + 1
    }
    if(itertest==200) stop("too many degenerate subsamples")
    
    Smat <- crossprod(res)
    Cmat <- det(Smat)^(-1/q)*Smat
    if(p>0)
      beta <- beta[2:(p+1),]
    else 
      beta <- matrix(0,0,q)
    
    if(k>0){
      tmp <- resmultiGS(X,Y,beta,Cmat,k,b,cc,0,convTol)
      gammarw <- tmp$gammarw
      scalerw <- tmp$scalerw
      betarw <- tmp$betarw
      rdisrw <- tmp$rdis
    }else{
      gammarw <- Cmat
    resrw <- res
    betarw <- beta
    places <- t(combn(1:n,2))
    term1vector <- as.matrix(res[places[,1],])
    term2vector <- as.matrix(res[places[,2],])
    diffrw <- term1vector-term2vector
    
    rdisrw <- sqrt(mahalanobis(diffrw, rep(0,q), gammarw))
    
    scalerw <- median(abs(rdisrw))/0.6745
    }
    
    if(i > 1){   
      scaletest <- mean(rhobiweight(rdisrw/sworst,cc))
      if(scaletest < b){ 
        
        ss <- sort(bestscales, index.return=TRUE)
        ind <- ss$ix[bestr]
        
        bestscales[ind] <- scalemultiGS(rdisrw,b,cc,scalerw)
        bestbetas[,ind] <- vecop(betarw)
        bestgammas[,ind] <- vecop(gammarw)
        sworst <- max(bestscales)
      }
    }else{
      bestscales[bestr] <- scalemultiGS(rdisrw,b,cc,scalerw)
      bestbetas[,bestr] <- vecop(betarw)
      bestgammas[,bestr] <- vecop(gammarw)
      
    }       
  }

  ibest <- which.min(bestscales)
  superbestscale <- bestscales[ibest]
  superbestbeta <- reconvec(bestbetas[,ibest],q)
  superbestgamma <- reconvec(bestgammas[,ibest],q)
  
  for(i in bestr:1){ 
    tmp <- resmultiGS(X,Y,reconvec(bestbetas[,i],q), reconvec(bestgammas[,i],q),maxIt,b,cc,bestscales[i],convTol)
    if(tmp$scalerw < superbestscale)  {
      superbestscale <- tmp$scalerw
      superbestgamma <- tmp$gammarw
      superbestbeta <- tmp$betarw
    }
  }
  GSbeta <- superbestbeta
  GSbeta <- as.matrix(GSbeta,drop=FALSE)
  
  return(list(betaR=t(GSbeta)))
}
