


##################################################################
##### Performance Comparison: Figure 8 and 11
##################################################################



# Use the following wrapper function to build models. 
# The input nsim corresponds to the index for snr.
# Each time wrapper is run, simulations with snr[nsim] 
# and parameters given in the wrapper function and are 
# repeated nrep times. 
# The same tuning procedure as in the Hastie et al. 2017
# is used.
# Parameters can be modified for different settings, with the 
# following for the low dimensional setting.

wrapper <- function(nsim){
  seed      <- 0
  
  # setting parameters
  n         <- 100      # training size
  p         <- 10       # feature dimension
  s         <- 5        # the number of signal features
  rho       <- 0.35     # autocorrelation between features
  beta.type <- 2        # patterns for beta
  snr       <- exp(seq(log(0.05),log(6),length=10))
  
  nval      <- 1000     # validation size
  nrep      <- 100      # simulation repetitions
  intercept <- FALSE    # including intercept or not in the model
  
  nlam      <- 50       # the number of lambda tuned in lasso and relaxed lasso
  nrelax    <- 10       # the number of gamma tuned in relaxed lasso
  maxsteps  <- 10       # the number of steps in fs and RandFS
  
  mtry <- ceiling(p*seq(1,10)/10) # mtry tuned in RandFS
  
  
  reg.funs  <- list()
  reg.funs[["Lasso"]]         <- function(x,y) lasso(x, y, intercept = intercept, nlam = nlam)
  reg.funs[["Stepwise"]]      <- function(x,y) fs(x, y, maxsteps = maxsteps, intercept = intercept)
  reg.funs[["Relaxed lasso"]] <- function(x,y) lasso(x, y, intercept = intercept, nrelax = nrelax, nlam = nlam)
  reg.funs[["RandFwd1"]]      <- function(x,y) randFwd(x, y, ntree = 1, maxsteps = maxsteps, mtry = p, intercept = intercept)
  reg.funs[["BagFwd"]]        <- function(x,y) randFwd(x, y, maxsteps = maxsteps,mtry=p, intercept = intercept)
  reg.funs[["RandFwd"]]       <- function(x,y) randFwd(x, y, maxsteps = maxsteps,intercept = intercept)
  reg.funs[["RandFwd mtry"]]  <- function(x,y) randFwd.mtry(x, y, maxsteps = maxsteps,intercept = intercept,mtry = mtry)
  
  sim.obj <- sim.master(n = n, p = p, nval = nval, reg.funs = reg.funs, nrep = nrep,
                        seed = seed, verbose = TRUE, rho = rho, s = s, beta.type = beta.type, 
                        file = NULL, file.rep = 0, snr = snr[nsim])
  
  if (nsim == 1) {
    result <- c(sim.obj,n = n, p = p, nval = nval, nrep = nrep, intercept = intercept,
                nlam = nlam, nrelax = nrelax, maxsteps = maxsteps, reg.funs = reg.funs)
  }else{
    result <- sim.obj
  }
  
  return(result)
}

# Apply wrapper to different snr
sim.obj     <- lapply(1:10, wrapper)
method.nums <- c(1,2,3,5,6,7)

# Generate the plot.
gp <- plot.result(sim.obj, snr, method.nums = method.nums, what = "error", type = "ave", std = T, distinguish = "point")
gp <- g1+theme(legend.just = c(0,1), legend.position = c(0.01,0.99))+
  ggtitle(paste(set.name[i], "Setting"))+
  scale_color_discrete(labels = c("BaggFS","Lasso","RandFS: mtry=0.33","RandFS","Relaxed lasso","FS"))+
  scale_shape_discrete(labels = c("BaggFS","Lasso","RandFS: mtry=0.33","RandFS","Relaxed lasso","FS"))+
  theme(plot.title       = element_text(hjust = 0.5,size=rel(2.5)), 
        legend.key.size  = unit(0.06, "npc"), 
        legend.spacing.y = unit(.02, "npc"),
        legend.text      = element_text(size=rel(2)),
        legend.title     = element_blank(),
        axis.title       = element_text(hjust=0.5,size=rel(2)),
        axis.text        = element_text(size=rel(1.7)))


##################################################################
##### Degrees of Freedom: Figure 7
##################################################################

# Simulations for estimating dof are mainly based on the R file 
# fig.df.R provided by 
# https://github.com/ryantibs/best-subset/blob/master/sims/fig.df.R



wrapper_df <- function(r){
  
  cat(r,"... ")
  eps = rnorm(n)*sigma
  y   = mu + eps
  
  # Lasso and Relaxed Lasso
  beta.las = coef(lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel))
  nzs.las  = colSums(beta.las != 0)[0:(nlam-1)*nrel+1]
  j = nlam - rev(match(p:0, rev(nzs.las), NA))
  ind = rep(j*nrel, each=nrel) + rep(1:nrel, p+1)
  yhat.las = (x %*% beta.las)[,ind]
  ip.las = colSums(yhat.las * eps)
  
  # Forward Selection
  yhat.fs = predict(fs(x,y,intercept=FALSE))
  ip.fs = colSums(yhat.fs * eps)[-(p+2)]
  
  # BaggFS
  yhat.bagfs <-  predict(randFwd(x,y,maxsteps=maxsteps,mtry=p,intercept = intercept))
  ip.bagfs <-  colSums(yhat.bagfs*eps)[-(p+2)]
  
  # RandFS with default mtry = p/3
  yhat.randfs <- predict(randFwd(x,y,maxsteps=maxsteps,intercept = intercept))
  ip.randfs <- colSums(yhat.randfs*eps)[-(p+2)]
  
  # RandFS with mtry = 0.25*p
  yhat.randfs_25 <- predict(randFwd(x,y,maxsteps=maxsteps,mtry=0.25*p,intercept = intercept))
  ip.randfs_25 <- colSums(yhat.randfs_25*eps)[-(p+2)]
  
  # RandFS with mtry = 0.1*p
  yhat.randfs_10 <- predict(randFwd(x,y,maxsteps=maxsteps,mtry=0.1*p,intercept = intercept))
  ip.randfs_10 <- colSums(yhat.randfs_10*eps)[-(p+2)]
  
  # RandFS with mtry = 0.6*p
  yhat.randfs_60 <- predict(randFwd(x,y,maxsteps=maxsteps,mtry=0.6*p,intercept = intercept))
  ip.randfs_60 <- colSums(yhat.randfs_60*eps)[-(p+2)]
  
  
  sim.obj <- enlist(n,p,s,rho,beta.type,snr,nval,seed,intercept,nlam,nrel,maxsteps,mu,sigma,
                    yhat.las,ip.las,yhat.fs,ip.fs,yhat.bagfs,ip.bagfs,yhat.randfs,ip.randfs,
                    yhat.randfs_10,ip.randfs_10,yhat.randfs_25,ip.randfs_25,
                    yhat.randfs_60,ip.randfs_60,r)

  return(sim.obj)
}



# Generate the dataset for estimating dof
seed <- 1
n    <- 70
p    <- 30
s    <- 5
rho  <- 0.35
beta.type <- 2
snr  <- 0.7
nval <- n
intercept <- FALSE
nlam <- 300
nrel <- 9
maxsteps <- 30

set.seed(seed)
xy.obj = sim.xy(n,p,nval,rho=rho,s=s,beta.type=beta.type,snr=snr)
x      = xy.obj$x
y      = xy.obj$y
mu     = as.numeric(x %*% xy.obj$beta)
sigma  = xy.obj$sigma



# Run the simulation 500 times to estimate dof
nrep <- 500
FinalResult <- lapply(1:nrep,wrapper_df)



# Aggeragate Data 
ip.las   <- matrix(0,nrep,(p+1)*nrel)
ip.fs    <- matrix(0,nrep,p+1)
ip.bagfs <- ip.randfs <- ip.randfs_25 <- 
  ip.randfs_10 <- ip.randfs_60 <- matrix(0,nrep,p+1)

for (r in 1:nrep) {
  ip.las[r,]       <-  FinalResult[[r]]$ip.las
  ip.fs[r,]        <-  FinalResult[[r]]$ip.fs
  ip.bagfs[r,]     <-  FinalResult[[r]]$ip.bagfs
  ip.randfs[r,]    <- FinalResult[[r]]$ip.randfs
  ip.randfs_10[r,] <- FinalResult[[r]]$ip.randfs_10
  ip.randfs_25[r,] <- FinalResult[[r]]$ip.randfs_25
  ip.randfs_60[r,] <- FinalResult[[r]]$ip.randfs_60
}

df.las       <-  colMeans(ip.las, na.rm=TRUE) / sigma^2
df.las       <-  matrix(df.las, p+1, nrel, byrow=TRUE)
df.fs        <-  colMeans(ip.fs, na.rm=TRUE) / sigma^2
df.bagfs     <- colMeans(ip.bagfs,na.rm=TRUE)/sigma^2
df.randfs    <- colMeans(ip.randfs,na.rm=TRUE)/sigma^2
df.randfs_10 <- colMeans(ip.randfs_10,na.rm=TRUE)/sigma^2
df.randfs_25 <- colMeans(ip.randfs_25,na.rm=TRUE)/sigma^2
df.randfs_60 <- colMeans(ip.randfs_60,na.rm=TRUE)/sigma^2



# The Left Plot in Figure 7
dat.df_left = data.frame(x      = rep(0:p,6),
                         y      = c(df.fs, df.bagfs, df.randfs, 
                                    df.las[,1], df.las[,5],df.las[,9]), 
                         Method = factor(rep(c("FS", "BaggFS","RandFS: mtry=0.33",
                                         "Lasso",
                                         "Relaxed lasso: \u03b3=0.5","Relaxed lasso: \u03b3=0"
                                         ), rep(p+1,6)))
                         )


gp.df_left <- ggplot(dat.df_left, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd = 0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd = 1) + geom_point( aes(shape = Method), size = 3) +
  theme_bw() + theme(legend.just = c(1,0), legend.pos = c(0.99,0.01))+
  theme(plot.title       = element_text(hjust = 0.5,size = rel(2.5)), 
        legend.key.size  = unit(0.06, "npc"), #legend.key = element_rect(fill = "gray88"),
        legend.spacing.y = unit(.02, "npc"),
        legend.text      = element_text(size = rel(2)),
        legend.title     = element_blank(),
        axis.title       = element_text(hjust = 0.5, size = rel(2)),
        axis.text        = element_text(size = rel(1.7)))



# The Right Plot in Figure 7
dat.df_right = data.frame(x      = rep(0:p,5),
                          y      = c(df.bagfs, df.randfs, df.randfs_10, 
                                   df.randfs_25, df.randfs_60), 
                          Method = factor(rep(c("BaggFS", "RandFS: mtry=0.33", 
                                              "RandFS: mtry=0.10","RandFS: mtry=0.25",
                                              "RandFS: mtry=0.60"),
                                            rep(p+1,5)))
                          )

gp.df_right <- ggplot(dat.df_right, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(aes(shape = Method), size = 3 ) +
  theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.99,0.01))+
  scale_color_manual(values=c("red","slateblue","deepskyblue","greenyellow","orange"))+
  theme(plot.title       = element_text(hjust = 0.5,size=rel(2.5)), 
        legend.key.size  = unit(0.06, "npc"), #legend.key = element_rect(fill = "gray88"),
        legend.spacing.y = unit(.02, "npc"),
        legend.text      = element_text(size=rel(2)),
        legend.title     = element_blank(),
        axis.title       = element_text(hjust=0.5,size=rel(2)),
        axis.text        = element_text(size=rel(1.7)))





















