
# This file contains plot functions, which are based
# on functions in the best-subset package.


### Plot the performance comparison results
plot.result= function(x,snr, method.nums=1:length(x$err.test), method.names=NULL,
                      what=c("error","risk","pve","nzero"), type=c("ave","med"), std=TRUE,
                      distinguish = c("color", "line", "point"),
                      lwd=1, pch=19, main=NULL, legend=TRUE, make.pdf=FALSE,
                      fig.dir=".", file.name="sim", w=6, h=6){
  
  # x :lists generated from the function sim.master.
  #    each list corresponds to one value of snr.
  #    x and snr must have the same length
  # distinguish: ways of distinguishing different lines in the plot.
  
  if (length(x) != length(snr)){
    stop("The number of columns of sim.obj must be the same as the length of snr!")
  }
  
  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }
  
  what = match.arg(what)
  type = match.arg(type)
  distinguish = match.arg(distinguish)
  if (is.null(method.names)) method.names = names(x[[1]]$err.test[method.nums])
  ii = method.nums
  
  # Populate the metric to be plotted
  N.s <- length(snr)
  err.rel.ave <- vector(mode="list",length=N.s)
  err.rel.std <- vector(mode="list",length=N.s)
  
  
  for(i in 1:N.s){
    sim.obj1 <- x[[i]]
    N = length(sim.obj1$err.test)
    
    # Calculating the metrics for tuning, relative error and relative risk
    # Temporary result for each SNR
    err.rel <-  vector(mode="list",length=N)
    
    if (what =="error"){
      for(j in 1:N){
        err.rel[[j]] = sim.obj1$err.test[[j]] / sim.obj1$sigma^2
      }
    } else if(what == "risk"){
      for(j in 1:N){
        err.rel[[j]]=sim.obj1$risk[[j]]/sim.obj1$risk.null
      }
    } else if(what == "nzero"){
      err.rel <- sim.obj1$nzs
    } else{
      err.rel <- sim.obj1$prop
    }
    
    err.obj <- tune.and.aggregate(sim.obj1, err.rel)
    
    if (type =="ave"){
      err.rel.ave[[i]] <- err.obj$z.val.ave
      err.rel.std[[i]] <- err.obj$z.val.std
    } else{
      err.rel.ave[[i]] <- err.obj$z.val.med
      err.rel.std[[i]] <- err.obj$z.val.mad
    }
  }
  
  # modify the x labels 
  x_lab <- snr
  x_lab[seq(2, length(snr), by =2)] = ""
  
  # Produce the plot
  dat = data.frame(x=rep(snr,each=length(ii)),
                   y=as.numeric(matrix(unlist(err.rel.ave),nrow=N)[ii,]),
                   se=as.numeric(matrix(unlist(err.rel.std),nrow=N)[ii,]),
                   Method=factor(rep(method.names,N.s)))
  ylim <- range(min(dat$y)-dat$se[which.min(dat$y)],max(dat$y)+dat$se[which.max(dat$y)])
  # Generate ylab
  if (what =="error"){
    ylab0 <- "Relative test error (to Bayes)"
  } else if(what == "risk"){
    ylab0 <- "Relative risk (to null model)"
  } else if(what == "nzero"){
    ylab0 <- "Number of nonzeros"
  } else{
    ylab0 <- "Proportion of variance explained"
  }
  
  
  if (distinguish == "color") {
    gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
      xlab("Signal-to-Noise Ratio") +
      ylab(ylab0) +
      geom_line(lwd=lwd) + geom_point(pch=pch) + theme_bw()+
      scale_x_continuous(trans="log",breaks=snr,
                         labels=x_lab)
  } else if(distinguish == "line"){
    gp = ggplot(dat, aes(x=x,y=y,color=Method, linetype = Method)) +
      xlab("Signal-to-Noise Ratio") +
      ylab(ylab0) +
      geom_line(lwd=lwd) + geom_point(pch=pch) + theme_bw()+
      scale_x_continuous(trans="log",breaks=snr,
                         labels=x_lab)
  } else{
    gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
      xlab("Signal-to-Noise Ratio") +
      ylab(ylab0) +
      geom_line(lwd=lwd) + geom_point( aes(shape = Method), size =3 ) + theme_bw()+
      scale_x_continuous(trans="log",breaks=snr,
                         labels=x_lab)
  }
  
  
  
  
  
  
  if (what=="pve") gp = gp + geom_line(aes(x=x, y=x/(1+x)), lwd=0.5,
                                       linetype=3, color="black")
  # if (what=="error") gp = gp + geom_line(aes(x=x, y=1+x), lwd=0.5,
  #                                        linetype=3, color="black")
  if (what =="nzero") gp = gp + geom_line(aes(x=x, y=sim.obj[[1]]$s), lwd=0.5,
                                          linetype=3, color="black")
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.2)
  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!legend) gp = gp + theme(legend.pos="none")
  if (!is.null(ylim)) gp = gp + coord_cartesian(ylim=ylim)
  if (make.pdf) {ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                        height=h, width=w, device="pdf")
  }else {gp
  }
}



# Plots for relative risk or error v.s. number of nonzero coefficients
# x-axis for randfwd or bagfwd is adjusted such that it represents the number 
# of nonzero coefficients in the base forward
plot.sim.randfwd <-   function(x, method.nums=1:length(x$err.test), method.names=NULL,
                               what=c("error","risk"), type=c("ave","med"), std=TRUE,
                               lwd=1, pch=19, main=NULL, legend=TRUE, make.pdf=FALSE,
                               fig.dir=".", file.name="sim", w=6, h=6) {
  
  
  
  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }
  
  what = match.arg(what)
  type = match.arg(type)
  if (is.null(method.names)) method.names = names(x$err.test[method.nums])
  ii = method.nums
  
  # Construct relative test error or relative risk
  N = length(x$err.test) # Number of methods
  err.rel = vector(mode="list",length=N)
  for (j in 1:N) {
    if (what=="error") err.rel[[j]] = x$err.test[[j]]/x$sigma[1]^2
    else err.rel[[j]] = x$risk[[j]]/x$risk.null
  }
  err.obj = tune.and.aggregate(x, err.rel)
  nzs.obj = tune.and.aggregate(x, x$nzs)
  names(nzs.obj$z.ave) <- names(nzs.obj$z.med) <- names(x$err.test)
  
  # obtain the effective nzs for randFwd, which is the same as the one in randFwd1
  nzs.obj$z.ave[["BagFwd"]] <- nzs.obj$z.ave[["RandFwd"]]  <- c(seq(0,x$maxsteps,by=1),x$maxsteps)
  nzs.obj$z.med[["BagFwd"]] <- nzs.obj$z.med[["RandFwd"]] <- c(seq(0,x$maxsteps,by=1),x$maxsteps)
  
  if (type=="ave") {
    xlist = nzs.obj$z.ave
    ylist = err.obj$z.ave
    ybars = err.obj$z.std
  }
  else {
    xlist = nzs.obj$z.med
    ylist = err.obj$z.med
    ybars = err.obj$z.mad
  }
  
  # Produce the plot
  dat = data.frame(x=unlist(xlist[ii]),
                   y=unlist(ylist[ii]),
                   se=unlist(ybars[ii]),
                   Method=factor(rep(method.names,
                                     lapply(xlist[ii],length))))
  
  
  gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
    xlab("Number of nonzero coefficients") +
    ylab(ifelse(what=="error",
                "Relative test error (to Bayes)",
                "Relative risk (to null model)")) +
    geom_line(lwd=lwd) + geom_point(pch=pch) + theme_bw()
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.2)
  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!legend) gp = gp + theme(legend.pos="none")
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
}


# Extract legends for grid arragement of plots 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}