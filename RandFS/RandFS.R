
# This file contains functions for randomized forward selection (RandFS)


#################################################################
### Random Forward
#################################################################

# ranndFwd is based on the function fs in the package best-subset.

# x: matrix of features
# y: vector of response
# maxsteps: the max number of steps in each forward selection
# intercept: Whether an intercept should be included
# normalize: whether predictors should be normalized before computing the path
# verbose: whether to print intermediate progress
# ntrees: the number of base random forward selections to be built
# mtry: the number of candidates for selection at each step
# replace: bootstrapping(T) or subsampling(F)
# samplesize: the size of bootstrap samples or subsamples

randFwd = function(x, y, maxsteps=min(nrow(x)-intercept,ncol(x),2000),
                   intercept=TRUE, normalize=TRUE, verbose=FALSE,
                   ntree=500,mtry=ceiling(ncol(x)/3),replace=TRUE,
                   samplesize= ifelse(replace,nrow(x),ceiling(.632*nrow(x) ) )) {
  
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
  # Save original x and y
  x0 = x
  y0 = y
  
  # Populate empty arrays for outputs from each random Fwd
  
  beta.all <- action.all <- df.all <- beta.all <- completepath.all <- 
    bls.all <- bx.all <- by.all <- vector(mode = "list",ntree)
  
  count.badtr <- 0
  # Start to build randomized fwd
  
  tr <- 1
  while(tr<=ntree){
    
    try.test <- try({train <- sample(1:n,samplesize,replace=replace)
    x.t <- x[train,,drop=FALSE]
    y.t <- y[train,drop=FALSE]
    
    # Center and scale, etc, the bootstrap sample.
    obj = standardize(x.t,y.t,intercept,normalize)
    x.t = obj$x
    y.t = obj$y
    bx.t = obj$bx
    by.t = obj$by
    sx.t = obj$sx
    
    # Randomly select a set of candidate features
    id.try <- sample(1:p,min(mtry,p),replace=FALSE) 
    
    
    #####
    # Find the first variable to enter and its sign
    z.t = scale(x.t[,id.try,drop=FALSE],center=F,scale=sqrt(colSums(x.t[,id.try,drop=FALSE]^2)))
    u.t = t(z.t) %*% y.t
    j.hit.try <- which.max(abs(u.t))  # Hitting coordinate in id.try
    j.hit = id.try[j.hit.try]        # Hitting coordinate of the feature
    sign.hit = Sign(u.t[j.hit.try])   # Hitting sign
    
    if (verbose) {
      cat(sprintf("1. Added variable %i, |A|=%i...",j.hit,1))
    }
    
    
    # Now iterate to find the sequence of FS estimates
    
    # Things to keep track of, and return at the end
    buf = min(maxsteps+1,500)
    action = numeric(buf)      # Actions taken
    df = numeric(buf)          # Degrees of freedom
    beta = matrix(0,p,buf)     # FS estimates
    
    
    # Record action, df, solution (df and solution are here
    # correspond to step 0; always a step behind)
    action[1] = j.hit
    df[1] = 0
    beta[,1] = 0
    
    
    # Other things to keep track of, but not return
    r = 1                       # Size of active set
    A = j.hit                   # Active set
    I = Seq(1,p)[-j.hit]        # Inactive set
    sign = sign.hit             # Active signs
    X1.t = x.t[,j.hit,drop=FALSE]   # Matrix X[,A] based on the bootstrap sample
    X2.t = x.t[,-j.hit,drop=FALSE]  # Matrix X[,I] based on the bootstrap sample
    k = 2                       # Step counter
    
    # Compute a skinny QR decomposition of X1.t
    qr.obj = qr(X1.t)
    Q = qr.Q(qr.obj,complete=TRUE)
    Q1 = Q[,1,drop=FALSE];
    Q2 = Q[,-1,drop=FALSE]
    R = qr.R(qr.obj)
    
    # Throughout the algorithm, we will maintain
    # the decomposition X1 = Q1*R. Dimensions:
    # X1: n x r
    # Q1: n x r
    # Q2: n x (n-r)
    # R:  r x r
    
    while (k <= maxsteps) {
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(action)) {
        buf = length(action)
        action = c(action,numeric(buf))
        df = c(df,numeric(buf))
        beta = cbind(beta,matrix(0,p,buf))
      }
      
      # Key quantities for the next entry
      a = backsolve(R,t(Q1) %*% y.t)      # Coefficients for the last step
      b = backsolve(R,t(Q1) %*% X2.t)     
      X2.resid.t = X2.t - X1.t %*% b      # Residuals of regressing on the active set
      
      
      # Draw candidate features
      mtry.k <- min(mtry, ncol(X2.resid.t))      # Trim mtry if needed
      id.try <- sample(1:ncol(X2.resid.t), mtry.k, replace=FALSE)
      
      z.t = scale(X2.resid.t[,id.try,drop=FALSE], center=F, scale=sqrt(colSums(X2.resid.t[,id.try,drop=FALSE]^2)))
      u.t = as.numeric(t(z.t) %*% y.t)
      
      
      # Otherwise find the next hitting time
      sign.u = Sign(u.t)
      abs.u = sign.u * u.t
      j.hit.try = which.max(abs.u)
      j.hit = id.try[j.hit.try]
      sign.hit = sign.u[j.hit.try]
      
      
      # Record action, df, solution
      action[k] = I[j.hit]
      df[k] = r
      beta[A,k] = a
      
      # Update rest of the variables
      r = r+1
      A = c(A,I[j.hit])
      I = I[-j.hit]
      sign = c(sign,sign.hit)
      X1.t = cbind(X1.t,X2.t[,j.hit])
      X2.t = X2.t[,-j.hit,drop=FALSE]
      
      # Update the QR decomposition
      updated.qr = updateQR(Q1,Q2,R,X1.t[,r])
      Q1 = updated.qr$Q1
      Q2 = updated.qr$Q2
      R = updated.qr$R
      
      if (verbose) {
        cat(sprintf("\n%i. Added variable %i, |A|=%i...",k,A[r],r))
      }
      
      # Update counter
      k = k+1
      k
    }
    
    
    # Record df and solution at last step
    df[k] = k-1
    beta[A,k] = backsolve(R,t(Q1) %*% y.t)
    
    # Trim
    action = action[Seq(1,k-1)]
    df = df[Seq(1,k)]
    beta = beta[,Seq(1,k),drop=FALSE]
    
    # If we stopped short of the complete path, then note this
    # Else we computed the complete path, so record LS solution
    if (k-1 < min(n-intercept,p)) {
      completepath = FALSE
      bls = NULL
    }else {
      completepath = TRUE
      bls = beta[,k]
    }
    
    if (verbose) cat("\n")
    
    # Adjust for the effect of centering and scaling
    if (intercept) df = df+1
    if (normalize) beta = beta/sx.t
    if (normalize && completepath) bls = bls/sx.t
    
    # Assign column names
    colnames(beta) = as.character(Seq(0,k-1))
    
    action.all[[tr]] <- action
    df.all[[tr]] <- df
    beta.all[[tr]] <- beta
    completepath.all[[tr]] <- completepath
    bls.all[[tr]] <- bls
    bx.all[[tr]] <- bx.t
    by.all[[tr]] <- by.t
    
    if (verbose) cat(paste(tr, " random forward selection completed.\n"))
    },
    
    silent=T
    )
    if(class(try.test) != "try-error"){
      tr <- tr+1
    }
    
  }
  
  
  action <- Reduce("+",action.all)/length(action.all)
  df <- Reduce("+",df.all)/length(df.all)
  beta <- Reduce("+",beta.all)/length(beta.all)
  completepath <- Reduce("+",completepath.all)/length(completepath.all)
  completepath <- as.logical(completepath)
  bls <- Reduce("+",bls.all)/length(bls.all)
  
  # calculate center, scale, etc of original x and y 
  obj.temp = standardize(x0,y0,intercept,normalize)
  bx = obj.temp$bx
  by = obj.temp$by
  sx = obj.temp$sx
  
  
  
  #### action, df, bx, by here are not meanful
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize,
             action.all=action.all,df.all=df.all,beta.all=beta.all,
             completepath.all=completepath.all,bls.all=bls.all,bx.all=bx.all,by.all=by.all)
  
  class(out) = "randFwd"
  
  return(out)
}



###### randfwd with mtry as a tuning parameter  ####
# The following is modified based on randFwd above 
# so that both mtry and step k are considered as tuning parameters.



randFwd.mtry = function(x, y, maxsteps=min(nrow(x)-intercept,ncol(x),2000),
                        intercept=TRUE, normalize=TRUE, verbose=FALSE,
                        ntree=500,mtry=ceiling(ncol(x)/3),replace=TRUE,
                        samplesize= ifelse(replace,nrow(x),ceiling(.632*nrow(x) ) )) {
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
  # Save original x and y
  x0 = x
  y0 = y
  
  # Populate empty lists for storing results from each mtry
  n.mtry <- length(mtry)
  action.mtry <- df.mtry <- beta.mtry <- completepath.mtry <- bls.mtry <- 
    bx.mtry <- by.mtry <- vector(mode="list",n.mtry)
  
  # loop through mtry
  for (i.mtry in 1:n.mtry) {
    
    # Populate empty arrays for outputs from each random Fwd
    
    beta.all <- action.all <- df.all <-  completepath.all <- 
      bls.all <- bx.all <- by.all <- vector(mode = "list",ntree)
    
    count.badtr <- 0
    # Start to build randomized fwd
    
    tr <- 1
    while(tr<=ntree){
      try.test <- try({train <- sample(1:n,samplesize,replace=replace)
      x.t <- x[train,,drop=FALSE]
      y.t <- y[train,drop=FALSE]
      
      # Center and scale, etc, the bootstrap sample.
      obj = standardize(x.t,y.t,intercept,normalize)
      x.t = obj$x
      y.t = obj$y
      bx.t = obj$bx
      by.t = obj$by
      sx.t = obj$sx
      
      # Randomly select a set of candidate features
      id.try <- sample(1:p,min(mtry[i.mtry],p),replace=FALSE) 
      
      
      #####
      # Find the first variable to enter and its sign
      z.t = scale(x.t[,id.try,drop=FALSE],center=F,scale=sqrt(colSums(x.t[,id.try,drop=FALSE]^2)))
      u.t = t(z.t) %*% y.t
      j.hit.try <- which.max(abs(u.t))  # Hitting coordinate in id.try
      j.hit = id.try[j.hit.try]        # Hitting coordinate of the feature
      sign.hit = Sign(u.t[j.hit.try])   # Hitting sign
      
      if (verbose) {
        cat(sprintf("1. Added variable %i, |A|=%i...",j.hit,1))
      }
      
      
      # Now iterate to find the sequence of FS estimates
      
      # Things to keep track of, and return at the end
      buf = min(maxsteps+1,500)
      action = numeric(buf)      # Actions taken
      df = numeric(buf)          # Degrees of freedom
      beta = matrix(0,p,buf)     # FS estimates
      
      
      # Record action, df, solution (df and solution are here
      # correspond to step 0; always a step behind)
      action[1] = j.hit
      df[1] = 0
      beta[,1] = 0
      
      
      # Other things to keep track of, but not return
      r = 1                       # Size of active set
      A = j.hit                   # Active set
      I = Seq(1,p)[-j.hit]        # Inactive set
      sign = sign.hit             # Active signs
      X1.t = x.t[,j.hit,drop=FALSE]   # Matrix X[,A] based on the bootstrap sample
      X2.t = x.t[,-j.hit,drop=FALSE]  # Matrix X[,I] based on the bootstrap sample
      k = 2                       # Step counter
      
      # Compute a skinny QR decomposition of X1.t
      qr.obj = qr(X1.t)
      Q = qr.Q(qr.obj,complete=TRUE)
      Q1 = Q[,1,drop=FALSE];
      Q2 = Q[,-1,drop=FALSE]
      R = qr.R(qr.obj)
      
      # Throughout the algorithm, we will maintain
      # the decomposition X1 = Q1*R. Dimensions:
      # X1: n x r
      # Q1: n x r
      # Q2: n x (n-r)
      # R:  r x r
      
      while (k <= maxsteps) {
        ##########
        # Check if we've reached the end of the buffer
        if (k > length(action)) {
          buf = length(action)
          action = c(action,numeric(buf))
          df = c(df,numeric(buf))
          beta = cbind(beta,matrix(0,p,buf))
        }
        
        # Key quantities for the next entry
        a = backsolve(R,t(Q1) %*% y.t)      # Coefficients for the last step
        b = backsolve(R,t(Q1) %*% X2.t)     
        X2.resid.t = X2.t - X1.t %*% b      # Residuals of regressing on the active set
        
        
        # Draw candidate features
        mtry.k <- min(mtry[i.mtry], ncol(X2.resid.t))      # Trim mtry if needed
        id.try <- sample(1:ncol(X2.resid.t),mtry.k,replace=FALSE)
        
        z.t = scale(X2.resid.t[,id.try,drop=FALSE],center=F,scale=sqrt(colSums(X2.resid.t[,id.try,drop=FALSE]^2)))
        u.t = as.numeric(t(z.t) %*% y.t)
        
        
        # Otherwise find the next hitting time
        sign.u = Sign(u.t)
        abs.u = sign.u * u.t
        j.hit.try = which.max(abs.u)
        j.hit <- id.try[j.hit.try]
        sign.hit = sign.u[j.hit.try]
        
        
        # Record action, df, solution
        action[k] = I[j.hit]
        df[k] = r
        beta[A,k] = a
        
        # Update rest of the variables
        r = r+1
        A = c(A,I[j.hit])
        I = I[-j.hit]
        sign = c(sign,sign.hit)
        X1.t = cbind(X1.t,X2.t[,j.hit])
        X2.t = X2.t[,-j.hit,drop=FALSE]
        
        # Update the QR decomposition
        updated.qr = updateQR(Q1,Q2,R,X1.t[,r])
        Q1 = updated.qr$Q1
        Q2 = updated.qr$Q2
        R = updated.qr$R
        
        if (verbose) {
          cat(sprintf("\n%i. Added variable %i, |A|=%i...",k,A[r],r))
        }
        
        # Update counter
        k = k+1
        k
      }
      
      
      # Record df and solution at last step
      df[k] = k-1
      beta[A,k] = backsolve(R,t(Q1) %*% y.t)
      
      # Trim
      action = action[Seq(1,k-1)]
      df = df[Seq(1,k)]
      beta = beta[,Seq(1,k),drop=FALSE]
      
      # If we stopped short of the complete path, then note this
      # Else we computed the complete path, so record LS solution
      if (k-1 < min(n-intercept,p)) {
        completepath = FALSE
        bls = NULL
      }else {
        completepath = TRUE
        bls = beta[,k]
      }
      
      if (verbose) cat("\n")
      
      # Adjust for the effect of centering and scaling
      if (intercept) df = df+1
      if (normalize) beta = beta/sx.t
      if (normalize && completepath) bls = bls/sx.t
      
      # Assign column names
      colnames(beta) = as.character(Seq(0,k-1))
      
      action.all[[tr]] <- action
      df.all[[tr]] <- df
      beta.all[[tr]] <- beta
      completepath.all[[tr]] <- completepath
      bls.all[[tr]] <- bls
      bx.all[[tr]] <- bx.t
      by.all[[tr]] <- by.t
      
      if (verbose) cat(paste(tr, " random forward selection completed.\n"))
      },
      
      silent=T
      )
      if(class(try.test) != "try-error"){
        tr <- tr+1
      }
    }
    
    
    action.mtry[[i.mtry]] <- Reduce("+",action.all)/length(action.all)
    df.mtry[[i.mtry]] <- Reduce("+",df.all)/length(df.all)
    beta.mtry[[i.mtry]] <- Reduce("+",beta.all)/length(beta.all)
    completepath.mtry[[i.mtry]] <- Reduce("+",completepath.all)/length(completepath.all)
    completepath.mtry[[i.mtry]] <- as.logical(completepath)
    bls.mtry[[i.mtry]] <- Reduce("+",bls.all)/length(bls.all)
    bx.mtry[[i.mtry]] <- Reduce("+",bx.all)/length(bx.all)
    by.mtry[[i.mtry]] <- Reduce("+",by.all)/length(by.all)
    
  }
  
  action.mtry.bind <- do.call(cbind,action.mtry)
  df.mtry.bind <- do.call(cbind,df.mtry)
  beta.mtry.bind <- do.call(cbind,beta.mtry)
  completepath.mtry.bind <- do.call(cbind,completepath.mtry)
  bls.mtry.bind <- do.call(cbind,bls.mtry)
  bx.mtry.bind <- do.call(cbind,bx.mtry)
  by.mtry.bind <- do.call(cbind,by.mtry)
  
  
  # calculate center, scale, etc of original x and y 
  obj.temp = standardize(x0,y0,intercept,normalize)
  bx = obj.temp$bx
  by = obj.temp$by
  sx = obj.temp$sx
  
  
  #### action, df, bx, by here are not meanful
  out = list(action=action.mtry.bind,df=df.mtry.bind,beta=beta.mtry.bind
             ,completepath=completepath.mtry.bind,bls=bls.mtry.bind,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize,mtry=mtry
  )
  
  class(out) = "randFwd"
  
  return(out)
}



# Since we only consider integer steps in our simulations, the same 
# coefficient function and prediction function can be used here. 
# Detailed explanations of the parameters can be found in the file 
# best_subset.R.

coef.randFwd = function(object, s, ...) {
  beta = object$beta
  if (object$completepath) beta = cbind(beta,object$bls)
  k = ncol(beta)-1
  if (missing(s)) s = 0:k
  else if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  knots = 0:k
  decreasing = FALSE
  
  beta.mat = coef.interpolate(beta,s,knots,decreasing)
  if (object$intercept) return(rbind(rep(object$by,ncol(beta.mat)),beta.mat))
  else return(beta.mat)
}


predict.randFwd = function(object, newx, s, ...) {
  beta = coef.randFwd(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}





