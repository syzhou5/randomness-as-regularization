
# This file contains functions used in estimating the 
# degrees of freedom (DoF) of random forests (RFs).



# Simulating data from the MARS model, which contains 
# interaction terms.
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

# Simulating data from the MARSadd model.
sim.mars2 <- function(n,nval,sigma){
  x <- matrix(runif(10*n,0,1),n,10)
  xval <- matrix(runif(10*nval,0,1),nval,10)
  e <- matrix(rnorm(n,0,sigma),n,1)
  eval <- matrix(rnorm(nval,0,sigma),nval,1)
  mu <- as.vector(0.1*exp(4*x[,1])+4/(1+exp(-20*(x[,2]-0.5)))+
                    3*x[,3]+2*x[,4]+x[,5]        )
  muval <- as.vector(0.1*exp(4*xval[,1])+4/(1+exp(-20*(xval[,2]-0.5)))+
                       3*xval[,3]+2*xval[,4]+xval[,5]        )
  y <- as.vector(0.1*exp(4*x[,1])+4/(1+exp(-20*(x[,2]-0.5)))+
                   3*x[,3]+2*x[,4]+x[,5] + e          )
  yval <- as.vector(0.1*exp(4*xval[,1])+4/(1+exp(-20*(xval[,2]-0.5)))+
                      3*xval[,3]+2*xval[,4]+xval[,5] + eval         )
  
  enlist(x,y,xval,yval,mu,muval)
}



# Simulating data from the linear model, which is the same as 
# the one in the R package best-subset.
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


# obtaining test mse and predictions from RFs.
RF.hat <- function(sim.obj,mtry=NULL,maxnodes=NULL,ntree=500){
  x <- sim.obj$x
  y <- sim.obj$y
  
  if (is.null(mtry)) {
    obj <- randomForest(x,y,xtest=x,ytest=y,maxnodes = maxnodes)
  }else{
    obj <- randomForest(x,y,xtest=x,ytest=y,mtry=mtry,maxnodes = maxnodes)
  }
  
  mse <- obj$test$mse[500]
  yhat <- obj$test$predicted
  return(enlist(mse,yhat))
}















