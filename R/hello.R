# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

gendata <- function(N, TT, seed, set = 1, r= 0.5){
  p <- 2
  vart=c(3,6,8,11,17,23)
  g1=1
  g2=2
  g3=3
  g4=4
  g5=2
  g6=3.5
  g7=2.5
  f=function(t)(g1*(t<vart[1])+g2*(t>=vart[1])*(t<vart[2])+g3*(t>=vart[2])*(t<vart[3])+g4*(t>=vart[3])*(t<vart[4])+g5*(t>=vart[4])*(t<vart[5])+g6*(t>=vart[5])*(t<vart[6])+g7*(t>=vart[6]))
  Time =seq(1,TT,length.out=TT)
  Ft <- f(Time)
  betat <- matrix(f(seq(TT+1, TT+p*TT)), nrow=TT, ncol=p)
  Alpha <- rnorm(N)
  epsi <- matrix(0, N, TT)
  Xit <- array(0, dim=c(N, TT, 2))
  set.seed(seed)
  if(set==1){
    for(i in 1:N){
      epsi[i, ] <- as.vector(arima.sim(list(order=c(1,0,0),ar=0.5),n=TT))
      Xit[i,,1] <- as.vector(arima.sim(list(order=c(1,0,0),ar=0.5),n=TT))
      Xit[i,,2] <- as.vector(arima.sim(list(order=c(1,0,0),ar=0.5),n=TT))
    }
  }else if(set==2){
    for(i in 1:N){
      epsi[i, ] <- as.vector(arima.sim(list(order=c(1,0,0),ar=0.5),n=TT))
      Xit[i,,1] <- as.vector(arima.sim(list(order=c(1,0,0),ar=0.5),n=TT))
      Xit[i,,2] <- r * Xit[i,,1] + sqrt(1-r^2) * rnorm(TT)
    }
  }

  Yit <- matrix(Ft, nrow=N, ncol=TT, byrow=T) +
    matrix(betat[,1], nrow=N, ncol=TT, byrow=T) * Xit[,,1] +
    matrix(betat[,2], nrow=N, ncol=TT, byrow=T) * Xit[,,2] +
    matrix(Alpha, nrow=N, ncol=TT) + epsi
  gt <- gamma1t <- gamma2t <- rep(0,TT)
  Ft <- Ft - Ft[1]
  gt[1]=Ft[1]
  gt[2:TT]=Ft[2:TT]-Ft[1:(TT-1)]
  gamma1t[1] <- betat[1,1]; gamma1t[2:TT] <- betat[2:TT, 1] - betat[1:(TT-1), 1]
  gamma2t[1] <- betat[1,2]; gamma2t[2:TT] <- betat[2:TT, 2] - betat[1:(TT-1), 2]
  return(list(Yit=Yit, Xit = Xit, gt0=gt, gamma1t0=gamma1t, gamma2t0 = gamma2t,
              ft0=Ft, beta1t0= betat[,1], beta2t0=betat[,2]))
}

constructXrow <- function(X, i, s, t){
  p <- dim(X)[3]
  TT <- dim(X)[2]
  z <- numeric(p*TT)
  for(j in 1:p){
    for( is in 1:t){
      if(is <= s){
        z[(j-1)*TT +is] <- X[i,t,j] - X[i,s,j]
      }else{
        z[(j-1)*TT +is] <- X[i,t,j]
      }
    }
  }
  return(z)
}

# nptvcm is simplification of NonParametric Time-Varying Coefﬁcients Model
nptvcm <- function(Xit, Yit, lambda=NULL, tuning.select='CV', nfolds=10, estimateSE=F, penalty='SCAD',
                   nlambda=50, lambda.max=5){
  require(ncvreg)
  N <- nrow(Yit); TT <- ncol(Yit)
  p <- dim(Xit)[3]
  X <- matrix(0, N*TT*(TT-1)/2, (p+1)*TT - 1)
  Tend <- ncol(X)
  k <- 1
  for( i in 1:N){
    for(s in 1:(TT-1)){
      #for(s in 1){
      for(t in (s+1):TT){
        X[k, s:(t-1)] =1
        X[k, TT:Tend] <- constructXrow(Xit, i, s, t)
        k <- k+1
      }
    }
  }
  Y <- numeric(N*(TT-1)/2)
  k <- 1
  for(i in 1:N){
    for(s in 1:(TT-1)){
      # for(s in 1){
      for(t in (s+1):TT){
        Y[k] <- Yit[i,t]-Yit[i,s]
        k <- k + 1
      }
    }
  }
  lambda.info <- list()
  lambda.info$tuning.select <- tuning.select
  if(is.null(lambda) && tuning.select=='CV'){
    cvfit <- cv.ncvreg(X, Y, nfolds = nfolds, penalty=penalty, nlambda=nlambda)
    lambda <- cvfit$lambda.min
    lambda.info$lambda.set  <- cvfit$lambda
    lambda.info$cv.error <- cvfit$cve
    lambda.info$lambda.minid <- cvfit$min
    lambda.info$select.lambda <- lambda
  }else if(is.null(lambda) && tuning.select=='BIC'){
    #nlambda <- 50; lambda.max <- 4
    lambda_set <- exp(seq(log(lambda.max), log(0.001 * lambda.max),
                          len = nlambda - 1))
    BICfunc <- function(Stheta, theta, N){
      log(Stheta) + 1/2*sum(theta!=0)*log(N)/N
    }
    nlam <- length(lambda_set)
    BICvalue <- numeric(nlam)
    for(ilam in 1:nlam){
      fit <- ncvreg(X, Y, lambda=lambda_set[ilam])
      epsi_pena <- Y - cbind(1,X) %*% fit$beta
      Stheta <- sum(epsi_pena^2)/N
      bic <- BICfunc(Stheta, fit$beta, N)
      BICvalue[ilam] <- bic
    }
    minid <- which.min(BICvalue)
    lambda <- lambda_set[minid]
    lambda.info$lambda.set <- lambda_set
    lambda.info$BICvalue <- BICvalue
    lambda.info$lambda.minid <- minid
    lambda.info$select.lambda <- lambda
  }
  fit <- ncvreg(X, Y, lambda=lambda, penalty=penalty)
  htheta_pena <- fit$beta[2:(ncol(X)+1)]
  res <- list()
  Tend <- (p+1)*TT -1
  W <- matrix(0, Tend, Tend)
  for(j in 1:(p+1)){
    if(j == 1){
      for(t in 1:(TT-1)){
        W[1:t, t] <- rep(1, t)
      }
    }else{
      for(t in 1:TT){
        W[((j-1)*TT): ((j-1)*TT+t-1), (j-1)*TT+t-1] <- rep(1, t)
      }
    }
  }
  betaf_pena <- c(0, t(W) %*% htheta_pena)
  betaf_pena_mat <- matrix(betaf_pena, TT, p+1)
  res$hbetaf <- betaf_pena_mat
  if(estimateSE){
    epsi_pena <- Y - cbind(1,X) %*% fit$beta
    index_pena <- which(htheta_pena!= 0)
    hLambdainv_pena <- matrix(0, Tend , Tend )
    hLambdainv_pena[index_pena, index_pena] <- qr.solve(1/N* t(X[,index_pena])%*% X[,index_pena])
    Weps_pena <- X * matrix(epsi_pena, nrow(X), ncol(X))
    Xarray <- array(0, dim=c(N, (p+1)*TT-1, TT*(TT-1)/2))
    for(i in 1:N){
      Xarray[i, ,] <- t(Weps_pena[((i-1)*TT*(TT-1)/2+1):(i*TT*(TT-1)/2),])
    }
    SWeps_pena <- apply(Xarray, c(1,2), sum)
    hGamma_pena <- matrix(0, Tend, Tend)

    hGamma_pena[index_pena, index_pena] <- 1/N * t(SWeps_pena[,index_pena]) %*% SWeps_pena[,index_pena]
    hSigma_theta_pena <- 1/N * hLambdainv_pena%*% hGamma_pena %*% hLambdainv_pena
    Var_raw_pena <- diag(t(W) %*% hSigma_theta_pena %*% W)
    se_raw_pena <- sqrt(Var_raw_pena)

    se_raw_pena_mat <- matrix(c(0, se_raw_pena), TT, p+1)
    res$se <- se_raw_pena_mat
  }
  res$selectlambda.info <- lambda.info
  class(res) <- 'nptvcm'
  return(res)
}
plot <- function(x, ...) UseMethod("plot")
plot.nptvcm <- function(res){
  if(!inherits(res, 'nptvcm')) stop('The class of res is wrong!')
  numfun <- ncol(res$hbetaf)
  TT <- nrow(res$hbetaf)
  par(mfrow=c(ceiling(numfun/2),2))
  if(res$selectlambda.info$tuning.select=='CV'){
    minid <- res$selectlambda.info$lambda.minid
    interval <- max(1,(minid-10)):min(minid+10, length(res$selectlambda.info$lambda.set))
    plot(res$selectlambda.info$lambda.set[interval], res$selectlambda.info$cv.error[interval], type='o',
         xlab=expression(lambda), ylab='CV.error', main=paste0('min.lambda=',  round(res$selectlambda.info$select.lambda,digits = 3)) )
  }else if(res$selectlambda.info$tuning.select=='BIC'){
    minid <- res$selectlambda.info$lambda.minid
    interval <- max(1,(minid-10)):min(minid+10, length(res$selectlambda.info$lambda.set))
    plot(res$selectlambda.info$lambda.set[interval], res$selectlambda.info$BICvalue[interval], type='o',
         xlab=expression(lambda), ylab='BIC',main=paste0('min.lambda=', round(res$selectlambda.info$select.lambda,digits = 3) ))
  }
  for(j in 1:numfun){
    funjt <- res$hbetaf[,j]
    sdfunjt <- res$se[,j]
    lbetf <- funjt - 1.96*sdfunjt
    ubetf <- funjt + 1.96*sdfunjt
    tmp <- paste0('beta_',j-1)
    ylabel <- ifelse(j==1, expression('f(t)'), tmp)
    plot(1:TT,funjt,xlim=c(1,TT),xlab="t", ylab=ylabel,type="l",lty=1,col="black",lwd=2,
         ylim=c(min(lbetf-0.5), max(ubetf)+0.5))
    lines(1:TT,lbetf,lty=3,col="gray",lwd=2)
    lines(1:TT, ubetf ,lty=3,col="gray",lwd=2)
    polygon(c(1:TT,TT:1),c(lbetf,rev(ubetf)),col="gray",density = 20, angle = -45,fillOddEven = TRUE)
    abline(h =0,lty=4, col='green')
  }

}
