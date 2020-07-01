

# This file contains examples for comparing the performance of RFs 
# with diffrent mtry and find the optimal mtry at each SNR using 
# simulated data.

# Performance comparison of default RF and bagging is just comparison
# between RF with mtry = 1/3 and mtry = 1, so the codes can be applied
# by building RFs with only these two values of mtry. 
# Here to be consistent with the paper, we denote mtry as a proportion 
# of features available at each split.
# In R package randomForest used below, mtry is the number of 
# candidate features.

# Source functions in the R file Performance_Comparison_functions.R.
# Load the randomForest, ggplot2 and tidyverse package.
library(randomForest)
library(ggplot2)
library(tidyverse)


# Variance of the signal in MARS model is given here.
v.mars  <- 50.82657  # MARS model



# The followings are two functions for repeating the simulation 
# nrep times under the linear model setting and the MARS setting 
# respectively.
# Each time we run these functions, 10 SNR levels are under consideration. 
# At each SNR, we generate data from the corresponding model with 
# specified parameters and obtain the test MSE of RFs built with 
# various values of mtry.

wrapper_linear <- function(nsim){
  
  # initialize parameters
  set.seed(nsim)
  
  n    <- 500       # training size
  nval <- 500       # test size
  p    <- 20        # feature dimension
  rho  <- 0.35      # autocorrelation between features 
  s    <- 10        # the number of signal features
  beta.type <- 2    # pattern of beta
  
  snr  <- round(exp(seq(log(0.05),log(6),length=10)),digits = 2)
  
  
  reg.fun <- list();
  for (i in 1:p) {
    reg.fun[[paste("Random Forest_",i,sep="")]] <- function(sim.obj)RF(sim.obj,mtry=i)
  }
  
  N <- length(reg.fun)
  
  # Initialize the testing error as a list of the length N, the number of methods.
  # Each element in err.test is a matrix of 1xlength(snr), 
  # storing the value for corresponding snr.
  err.test <- vector(mode="list",length = N)
  reg.names <- names(reg.fun)
  names(err.test) <- reg.names
  
  for (k in 1:N) {
    err.test[[k]] <- matrix(NA,nrow=1,ncol=length(snr))
  }
  
  
  for (j in 1:length(snr)) {
    cat('j = ', j,'\n')
    xy.obj <- sim.xy(n,p,nval,rho,s,beta.type,snr[j])
    
    
    for (i in 1:N) {
      cat('i=', i,"\n")
      err.test[[i]][1,j] <- reg.fun[[i]](xy.obj)
    }
    
  }
  
  if(nsim == 1){
    result <- list(n=n,p=p,nval=nval,rho=rho,s=s,beta.type=beta.type,snr=snr,reg.fun=reg.fun,err.test=err.test )
  }else{
    result <- list(err.test=err.test )
  }
  
  return(result)
}


wrapper_mars <- function(nsim){
  # initialize parameters
  set.seed(nsim)
  
  n    <- 500    # training size
  nval <- 500    # test size
  p    <- 5      # feature dimension
  
  snr  <- round(exp(seq(log(0.05),log(6),length=10)),digits = 2)
  
  reg.fun <- list();
  
  reg.fun[["Random Forest_1"]] <- function(sim.obj) RF(sim.obj,mtry=1)
  reg.fun[["Random Forest_2"]] <- function(sim.obj) RF(sim.obj,mtry=2)
  reg.fun[["Random Forest_3"]] <- function(sim.obj) RF(sim.obj,mtry=3)
  reg.fun[["Random Forest_4"]] <- function(sim.obj) RF(sim.obj,mtry=4)
  reg.fun[["Random Forest_5"]] <- function(sim.obj) RF(sim.obj,mtry=5)
  
  N <- length(reg.fun)
  
  # Initialize the testing error as a list of the length N, the number of methods.
  # Each element in err.test is a matrix of 1xlength(snr), 
  # storing the value for corresponding snr.
  err.test <- vector(mode="list",length = N)
  reg.names <- names(reg.fun)
  names(err.test) <- reg.names
  
  for (k in 1:N) {
    err.test[[k]] <- matrix(NA,nrow=1,ncol=length(snr))
  }
  
  
  for (j in 1:length(snr)) {
    
    sigma <- (v.mars/snr[j])^0.5
    xy.obj <- sim.mars1(n,nval,sigma)
    
    for (i in 1:N) {
      err.test[[i]][1,j] <- reg.fun[[i]](xy.obj)
    }
  }
  
  if(nsim == 1){
    result <- list(n=n,p=p,nval=nval,snr=snr,reg.fun=reg.fun,err.test=err.test )
  }else{
    result <- list(err.test=err.test )
  }
  
  return(result)
}


# The following is an example for the linear model setting.
# Simulations with the MARS model can be carried out similarly
# using wrapper_mars.

# Run wrapper_linear 500 times.
FinalResult <- lapply(1:500, wrapper_linear)


# Extract parameters specified in the wrapper function.
snr       <- FinalResult[[1]]$snr    
p         <- FinalResult[[1]]$p
lsnr      <- length(snr)                    
nrep      <- length(FinalResult)
reg.names <- names(FinalResult[[1]]$reg.fun)
N         <- length(reg.names)  # the number of RFs built


#### There are two ways of determining the optimal mtry.
# The following is the one in Figure 5 in the paper.
# Optimal mtry is the one corresponding to the minimum
# average test MSE across 500 simulations.



# err.test is a list of length N corresponding to the N RFs 
# models built.
# Each element in err.test  is a nrep x lsnr matrix with 
# (i,j) entry corresponding to the test MSE of at i-th 
# simulation under j-th SNR as above.
# err.ave is a lsnr x N matrix with (i,j) entry corresponding 
# to the average test MSE of j-th RF model across nrep 
# simulations under i-th SNR level.

err.test <- vector(mode="list",length=N)
err.ave <- matrix(0,lsnr,N)

for (i in 1:N) {
  err.test[[i]] <- matrix(0,nrep,lsnr)
  for (j in 1:nrep) {
    err.test[[i]][j,] <- FinalResult[[j]]$err.test[[i]]
  }
  err.ave[,i] <- colMeans(err.test[[i]])
}


opt.mtry <- apply(err.ave,1,which.min)

dat1 <- data.frame(y=opt.mtry,x=snr,method=factor(rep("random forest",lsnr) ))
x_lab <- snr
x_lab[seq(2,length(snr),by=2)] = ""

gp1 <- ggplot(data=dat1,aes(x=x,y=y))+
  xlab("Signal-Noise Ratio")+
  ylab("Optimal mtry")+
  ggtitle("Optimal mtry under Linear Model")+
  geom_line(lwd=1)+geom_point(pch=19)+
  scale_x_continuous(trans = "log",breaks=snr,labels = x_lab)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2.5)), 
        axis.title = element_text(hjust=0.5,size=rel(2.2)),
        axis.text  = element_text(size=rel(2)))




# The following is the second way of determining optimal mtry
# which corresponds to Figure 10 in the paper.
# First find the optimal mtry in each of the 500 simulations.
# The optimal mtry reported is the average of the 500 observations.

# opt.mtry is a nrep x lsnr matrix, storing the values of optimal mtry for each simulation 
opt.mtry_2 <- matrix(0,nrep,lsnr)

for (i in 1:nrep) {
  err.test_2 <- FinalResult[[i]]$err.test %>% unlist() %>% matrix(nrow=10) %>% t()
  opt.mtry_2[i,] <- apply(err.test_2,2,which.min)/p
}

opt.mtry.ave <- apply(opt.mtry_2,2,mean,na.rm=TRUE)
opt.mtry.std <- apply(opt.mtry_2,2,sd,na.rm=TRUE)/sqrt(colSums(!is.na(opt.mtry_2)))

dat2 <- data.frame(y=opt.mtry.ave %>% unlist() %>% as.numeric(),
                  x=snr,
                  sd=opt.mtry.std %>% unlist() %>% as.numeric())

gp2 <- ggplot(data=dat2,aes(x=x,y=y))+
  xlab("Signal-Noise Ratio")+
  ylab("Optimal mtry")+
  ggtitle("Linear Model")+
  geom_line(lwd=1)+geom_point(pch=19 )+
  geom_errorbar(aes(ymin=y-sd,ymax=y+sd,width=0.2))+
  scale_x_continuous(trans = "log",breaks=snr,labels = x_lab)+
  theme_bw()+theme(legend.position = c(0.83,0.25))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2.5)), 
        axis.title = element_text(hjust=0.5,size=rel(2.2)),
        axis.text = element_text(size=rel(2)))

gp2  





















