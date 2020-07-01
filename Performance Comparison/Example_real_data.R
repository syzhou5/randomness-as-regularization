

# This file contains examples for comparing the performance of RF 
# with default mtry = 1/p and bagging (which is RF with mtry = 1)
# on real data.

# The RealData.RData file contains the real data used in the following
# simulations. Observations with missing values have been removed.
# The last column of each dataframe is the response.

# Source functions in the R file Performance_Comparison_functions.R.
# Load the randomForest, ggplot2 and tidyverse package.

library(randomForest)
library(ggplot2)
library(tidyverse)


# The following is a wrapper function for comparing the performance
# of RF and bagging on real datasets. 
# Given a dataset, each time we run wrapper, we obtain the K-fold 
# CV error of bagging and RF with various levels of noise injected
# into the response.
# This process is repeated 500 times.


wrapper <- function(nsim, dat){
  
  n <- dim(dat)[1]       # data size
  p <- dim(dat)[2] -1    # feature dimension
  K <- 10                # the number of folds in CV
  
  # ptg is the proportion alpha in page 15 of the paper.
  # Variance of the noise injected into the reponse is alpha 
  # porportion of the sample variance of the original response.
  ptg <- c(0,0.01,0.05,0.1,0.25,0.5)    
  sigma <- (var(dat[,p+1])*ptg)^0.5
  
  set.seed(nsim)
  cv.sigma <- lapply(sigma,noise_rf,dat=dat,K=K)
  
  if(nsim==1){
    result <- list(n=n,p=p,K=K,ptg=ptg,sigma=sigma,dat=dat,nsim=nsim,cv.sigma=cv.sigma)
  }else{
    result <- list(cv.sigma=cv.sigma)
  }
  
  return(result)
  
}


# The following is an example with the boston housing data.
# Simulations on other datasets can be carried out similarly.

# Run wrapper_linear 500 times.
FinalResult <- lapply(1:500, wrapper, dat=boston)


# Extract parameters in the simulation specified in wrapper.
ptg    <- FinalResult[[1]]$ptg
ls     <- length(ptg)
N      <- length(FinalResult)
p      <- FinalResult[[1]]$p
var.d  <- var(FinalResult[[1]]$dat[,p+1])
cv.bag <- cv.rf <-  matrix(NA, ls, N)  # K-fold CV error of bagging and RF

for (i in 1:N) {
  cv.temp    <- matrix(unlist(FinalResult[[i]]$cv.sigma),ncol=2,byrow = T)
  cv.bag[,i] <- cv.temp[,1]
  cv.rf[,i]  <- cv.temp[,2]
}

cv.RTE <- (cv.bag - cv.rf)/var.d*100     # relative test erorr (RTE)
cv.ave <- rowMeans(cv.RTE)               # average RTE from the 500 simulations
cv.sd  <- apply(cv.RTE, 1, sd, na.rm = T)/sqrt(rowSums(!is.na(cv.RTE)))

# Shift the RTE so that the plots based on the shifted RTE start from the origin.
RTE_shift <- cv.ave - cv.ave[1]

dat <- data.frame(y  = RTE_shift, 
                  x  = ptg,
                  sd = cv.sd)

gp <- ggplot(data = dat, aes(x = x,y = y)) +
  xlab("var(noise)/var(y)") +
  ylab("Shifted Relative Test Error(%)") +
  geom_line(lwd = 1) + 
  geom_point(pch=19) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.01) +
  theme_bw()+geom_hline(yintercept = 0,color="black",linetype="dashed")+
  theme(
    axis.title = element_text(hjust=0.5,size=rel(2.5)),
    axis.text = element_text(size=rel(2)))


















