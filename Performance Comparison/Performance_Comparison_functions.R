

# This file contains functions for comparing performance of
# Random Forests (RFs) with different mtry both on synthetic 
# data and real data.




# Generate data from the linear model. This function is from
# the R package best-subset.
sim.xy = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1) {
  
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  
  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }
  else Sigma = diag(1,p)
  
  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  enlist(x,y,xval,yval,Sigma,beta,sigma)
}




# Generate data from MARS model
# simga is the standard deviation of the error term. which is
# calculated based on a given SNR.
sim.mars1 <- function(n,nval,sigma){
  x <- matrix(runif(5*n,0,1),n,5)
  xval <- matrix(runif(5*nval,0,1),nval,5)
  e <- matrix(rnorm(n,0,sigma),n,1)
  eval <- matrix(rnorm(nval,0,sigma),nval,1)
  mu <- as.vector(10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.05)^2 + 10*x[,4]+5*x[,5])
  muval <- as.vector(10*sin(pi*xval[,1]*xval[,2])+20*(xval[,3]-0.05)^2 + 10*xval[,4]+5*xval[,5])
  y <- as.vector(10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.05)^2 + 10*x[,4]+5*x[,5]+e)
  yval <- as.vector(10*sin(pi*xval[,1]*xval[,2])+20*(xval[,3]-0.05)^2 + 10*xval[,4]+5*xval[,5]+eval)
  enlist(x,y,xval,yval,mu,muval)
}




# Build RFs and obtain test MSE
# sim.obj is a list generated from sim.xy or sim.mars1.
RF <- function(sim.obj,mtry=NULL,ntree=500){
  x <- sim.obj$x
  y <- sim.obj$y
  x.test <- sim.obj$xval
  y.test <- sim.obj$yval
  
  if (is.null(mtry)) {
    obj <- randomForest(x,y,xtest=x.test,ytest=y.test)
  } else{
    obj <- randomForest(x,y,xtest=x.test,ytest=y.test,mtry=mtry)
  }
  
  mse <- obj$test$mse[500]
  return(mse)
}




enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}





# Given a dataset and sd of noise injected into the
# response, calculate the K-fold CV error of RFs.
noise_rf <- function(dat,sigma,K){
  
  n <- dim(dat)[1]
  p <- dim(dat)[2]-1
  id <- sample(1:n,n,replace=F)
  
  # Inject noise into response and store it in the last coloumn.
  e <- rnorm(n,0,sigma)
  dat[,p+2] <- dat[,p+1]+e
  
  
  # Obtain the test error for each fold
  # 1st row for bagging and 2nd row for rf
  cv.error <- matrix(0,2,K)
  
  for (k in 1:K) {
    if (k !=K) {
      loc <- ((k-1)*floor(n/K)+1):(k*floor(n/K))
    }else{
      loc <- ((k-1)*floor(n/K)+1):n
    }
    test <- id[loc]
    
    m <- list(mtry=c(p,floor(p/3)))
    obj <- lapply(X=m[[1]],FUN=randomForest,x=dat[-test,1:p],y=dat[-test,p+2],
                  xtest=dat[test,1:p],ytest=dat[test,p+2],ntree=500)
    cv.error[1,k] <- obj[[1]][["test"]][["mse"]][500]
    cv.error[2,k] <- obj[[2]][["test"]][["mse"]][500]
    
  }
  
  cv <- rowMeans(cv.error)
  return(cv)
}
















