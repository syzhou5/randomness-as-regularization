
# This file contains examples of estimting Degrees of Freedom (DoF)
# of random forests (RF) under various settings.

# Source functions in the R file DoF_functions first.

library(randomForest)
library(ggplot2)


# Variance of signals of MARS and MARSadd models are given as following.
v.mars  <- 50.82657  # MARS model
v.mars2 <- 6.29620   # MARSadd Model


##### Estimating DoF of RF under the MARSadd Setting #####

# Generate Data
set.seed(1)
n       <- 1000   # training size
p       <- 10     # feature dimension
nval    <- 1000   # 

snr     <- 3.52
sigma   <- (v.mars2/snr)^0.5

xy.obj  <- sim.mars2(n,nval,sigma)
x       <- xy.obj$x
xval    <- xy.obj$xval
y       <- xy.obj$y
yval    <- xy.obj$yval

mu      <- xy.obj$mu

maxnd   <- ceiling(seq(2,n/5,length.out = 9))

# wrapper for repeating the simulations
wrapper <- function(r){
  reg.fun <- list();
  reg.fun[["Bagging"]] <- function(sim.obj,maxnodes) RF.hat(sim.obj,mtry=p,maxnodes = maxnodes)
  reg.fun[["Random Forest_2p/3"]] <- function(sim.obj,maxnodes) RF.hat(sim.obj,mtry=ceiling(2*p/3),maxnodes = maxnodes )
  reg.fun[["Random Forest_p/3"]] <- function(sim.obj,maxnodes) RF.hat(sim.obj,mtry=NULL,maxnodes = maxnodes)
  reg.fun[["Random Forest_p/10"]] <- function(sim.obj,maxnodes) RF.hat(sim.obj,mtry=ceiling(p/10),maxnodes = maxnodes)
  
  reg.names <- names(reg.fun)
  N <- length(reg.fun)
  
  ip <- vector(mode="list",length=N)
  names(ip) <- reg.names
  
  lmax <- length(maxnd)+1
  
  for (i in 1:N) {
    ip[[i]] <- matrix(0,1,lmax)
  }
  
  eps <- rnorm(n)*sigma
  y <- mu+eps
  
  xy.obj[["y"]] <- y
  for (i in 1:N) {
    for (j in 1:(lmax-1)) {
      yhat <- reg.fun[[i]](xy.obj,maxnd[[j]])[[2]]
      ip[[i]][1,j] <- sum(yhat*eps)
    }
    yhat <- reg.fun[[i]](xy.obj,maxnodes = NULL)[[2]]
    ip[[i]][1,lmax] <- sum(yhat*eps)
    
  }
  
  result <- list(reg.fun=reg.fun,N=N,ip=ip,lmax=lmax)
  return(result)
}



# Run the wrapper 500 times for estimating DoF
FinalResult <- lapply(1:500, wrapper)



# Aggregate data
N <- FinalResult[[1]]$N      # the number of models built
nrep <- length(FinalResult)
lmax <- length(maxnd)

ip.all <- vector(mode="list",N)
df.all <- matrix(0,lmax,N)

reg.names <- names(FinalResult[[1]]$reg.fun)

lgd.names <- c("Random Forest: mtry=1","Random Forest: mtry=0.67","Random Forest: mtry=0.33","Random Forest: mtry=0.1")

names(ip.all) <- reg.names
for (i in 1:N) {
  ip.all[[i]] <- matrix(0,nrep,lmax)
  for (j in 1:nrep) {
    ip.all[[i]][j,] <- FinalResult[[j]]$ip[[i]][-(lmax+1)]
  }
  df.all[,i] <- colMeans(ip.all[[i]],na.rm=TRUE)/sigma^2
  
}



# Draw the plot of DoF v.s. maxnodes
dat <- data.frame(x=rep(maxnd,N),y=matrix(df.all,ncol=1),method=factor(rep(lgd.names,each=lmax)))
gp<- ggplot(dat,aes(x=x,y=y,color=method))+
  ylab("Degree of Freedom")+
  xlab("maxnodes")+
  ggtitle("MARSadd Model")+
  geom_line(lwd=1) + geom_point(aes(shape = method), size =3 ) +
  geom_line(aes(y=x),linetype="dotted")+
  theme_bw() +
  theme(legend.justification = c(1,0),legend.position = c(0.99,0.01))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2.2)), 
        legend.key.size = unit(0.06, "npc"), 
        legend.spacing.y = unit(.02, "npc"),
        legend.text = element_text(size=rel(1.5)),
        legend.title = element_blank(),
        axis.title = element_text(hjust=0.5,size=rel(2)),
        axis.text = element_text(size=rel(1.7)))


##### Estimating DoF of RFs under the linear model Setting #####

# The DoF of RFs under the linear model setting is estimated in
# the same way as above except that data is generated with sim.xy
# instead of sim.mars2.

# Generate Data
n         <- 100      # training size
p         <- 10       # feature dimension
nval      <- 1000     # validation size
s         <- 5        # the number of signal features
rho       <- 0.35     # autocorrelation between features
beta.type <- 2        # patterns for beta
snr       <- 3.52     # signal-to-noise ratio

xy.obj    <- sim.xy(n, p, nval, rho, s, beta.type, snr)
x         <- xy.obj$x
xval      <- xy.obj$xval
y         <- xy.obj$y
yval      <- xy.obj$yval

mu        <- xy.obj$mu

maxnd     <- ceiling(seq(2,n/5,length.out = 9))

# The rest are the same as above in the MARSadd setting.


