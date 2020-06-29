








#' Predict function for fs object.
#'
#' Predict the response from a new set of predictor variables, using the
#'   coefficients from a particular step of the forward stepwise path.
#'
#' @param object The fs path object, as produced by the fs function.
#' @param newx Matrix of new predictor variables at which predictions should
#'   be made; if missing, the original (training) predictors are used.
#' @param s The step (or vector of steps) of the path at which coefficients
#'   should be computed. Can be fractional, in which case interpolation is
#'   performed. If missing, then the default is use all steps of the passed
#'   fs object.
#' @param ... Other arguments (currently not used).
#'
#' @details Note that at s = 1, there is one nonzero coefficient, at
#'   s = 2, there are two nonzero coefficients, etc. (This differs from the
#'   parametrization used in the \code{coef.fs} function in the R package
#'   \code{selectiveInference}, as the latter function delivers s-1 nonzero
#'   coefficients at step s, and was written to be consistent with the
#'   natural parametrization for the least angle regression path.)
#'
#' @export predict.fs
#' @export

predict.randFwd = function(object, newx, s, ...) {
  beta = coef.randFwd(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}




# To be reworked .....
# coef.randFwdMtry = function(object, s, ...) {
#   beta.org <- object$beta
#   
#   # Calculate the coefficients corresponding to each mtry with coef.randFwd
#   # Then aggregate the result back to a matrix
#   dim.temp <- dim(beta.org) 
#   nc <- dim.temp[2] # number of models considered in each randFwd.mtry
#   nr <- dim.temp[1] # number of features in the linear model
#   
#   n.mtry <- length(object$mtry)
#   K <- nc/n.mtry # number of steps for each mtry
#   
#   beta.all <- matrix(NA,nrow = nr,ncol = nc)
#   for (i in n.mtry) {
#     
#   }
#   
#   
#   
# }

#' Predict function for fs object.
#'
#' Predict the response from a new set of predictor variables, using the
#'   coefficients from a particular step of the forward stepwise path.
#'
#' @param object The fs path object, as produced by the fs function.
#' @param newx Matrix of new predictor variables at which predictions should
#'   be made; if missing, the original (training) predictors are used.
#' @param s The step (or vector of steps) of the path at which coefficients
#'   should be computed. Can be fractional, in which case interpolation is
#'   performed. If missing, then the default is use all steps of the passed
#'   fs object.
#' @param ... Other arguments (currently not used).
#'
#' @details Note that at s = 1, there is one nonzero coefficient, at
#'   s = 2, there are two nonzero coefficients, etc. (This differs from the
#'   parametrization used in the \code{coef.fs} function in the R package
#'   \code{selectiveInference}, as the latter function delivers s-1 nonzero
#'   coefficients at step s, and was written to be consistent with the
#'   natural parametrization for the least angle regression path.)
#'
#' @export predict.fs
#' @export

# predict.randFwd = function(object, newx, s, ...) {
#   beta = coef.randFwd(object,s)
#   if (missing(newx)) newx = object$x
#   else newx = matrix(newx,ncol=ncol(object$x))
#   
#   newx = scale(newx,object$bx,FALSE)
#   if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
#   return(newx %*% beta)
# }







randFwd_BO = function(x, y,xval,yval, maxsteps=min(nrow(x)-intercept,ncol(x),2000),
                      intercept=TRUE, normalize=TRUE, verbose=FALSE,
                      ntree=500,mtry=ceiling(p/3),replace=TRUE,
                      samplesize= ifelse(replace,nrow(x),ceiling(.632*nrow(x) ) )) {
  # randFwd_BO stands for randomized forward selection with each base optimized
  # x: matrix of features
  # y: vector of response
  # maxsteps: the max number of steps in each forward selection
  # intercept: Whether an intercept should be included
  # normalize: whether predictors should be normalized before computing the path
  # verbose: whether to print intermediate progress
  # ntrees: the number of random forward selections to be built
  # mtry: the number of candidates for selection at each step
  # replace: bootstrapping(T) or subsampling(F)
  # samplesize: the size of bootstrap or subsample
  # xval and yval are used to choose the optimal base learner.
  
  
  
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
  # out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
  #            x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  
  
  # beta.all <- action.all <- df.all <- beta.all <- completepath.all <- 
  #   bls.all <- bx.all <- by.all <- vector(mode = "list",ntree)
  beta.all <- df.all <- 
    bx.all <- by.all <- vector(mode = "list",ntree)
  bases <- vector(mode="list",ntree)
  
  
  # call randFwd function to build each base learner
  
  for (i in 1:ntree) {
    bases[[i]] <- randFwd(x=x,y=y,maxsteps=maxsteps,intercept=intercept,normalize=normalize,
                          verbose=verbose,ntree=1,mtry=mtry,replace=replace,
                          samplesize=samplesize)
    
    # calculate validation error to choose the best step
    muhat.val <- as.matrix(predict(bases[[i]],xval))
    err.val <- colMeans((muhat.val - yval)^2)
    opt.step <- which.min(err.val)
    beta.all[[i]] <- bases[[i]]$beta[,opt.step]
    bx.all[[i]] <- bases[[i]]$bx
    by.all[[i]] <- bases[[i]]$by
    
  }
  
  
  
  
  # Reduce("+", my.list) / length(my.list)
  # action <- Reduce("+",action.all)/length(action.all)
  # df <- Reduce("+",df.all)/length(df.all)
  beta <- Reduce("+",beta.all)/length(beta.all)
  # completepath <- Reduce("+",completepath.all)/length(completepath.all)
  # completepath <- as.logical(completepath)
  # bls <- Reduce("+",bls.all)/length(bls.all)
  bx <- Reduce("+",bx.all)/length(bx.all)
  by <- Reduce("+",by.all)/length(by.all)
  
  
  
  #### action, df, bx, by here are not meanful
  out = list(bases=bases,  beta=beta,x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize,
             beta.all=beta.all,bx.all=bx.all,by.all=by.all)
  
  class(out) = "randFwd_BO"
  
  return(out)
}

coef.randFwd_BO = function(object,  ...) {
  beta = object$beta
  # k = length(beta)
  # if (missing(s)) s = 0:k
  # else if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  # knots = 0:k
  # decreasing = FALSE
  # 
  # beta.mat = coef.interpolate(beta,s,knots,decreasing)
  # if (object$intercept) return(rbind(rep(object$by,ncol(beta.mat)),beta.mat))
  # else return(beta.mat)
  if (object$intercept) return(c(object$by,beta))
  else return(beta)
}

#' Predict function for fs object.
#'
#' Predict the response from a new set of predictor variables, using the
#'   coefficients from a particular step of the forward stepwise path.
#'
#' @param object The fs path object, as produced by the fs function.
#' @param newx Matrix of new predictor variables at which predictions should
#'   be made; if missing, the original (training) predictors are used.
#' @param s The step (or vector of steps) of the path at which coefficients
#'   should be computed. Can be fractional, in which case interpolation is
#'   performed. If missing, then the default is use all steps of the passed
#'   fs object.
#' @param ... Other arguments (currently not used).
#'
#' @details Note that at s = 1, there is one nonzero coefficient, at
#'   s = 2, there are two nonzero coefficients, etc. (This differs from the
#'   parametrization used in the \code{coef.fs} function in the R package
#'   \code{selectiveInference}, as the latter function delivers s-1 nonzero
#'   coefficients at step s, and was written to be consistent with the
#'   natural parametrization for the least angle regression path.)
#'
#' @export predict.fs
#' @export

predict.randFwd_BO = function(object, newx,  ...) {
  beta = coef.randFwd_BO(object)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}







sim.master2 = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                       file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1,indicator) {
  
  # indicator is a vector whose TRUE means reg is randFwd_BO
  # so indictor must has the same length as reg.funs
  # cat("checkpoint1:",p)
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  
  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)
  
  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = 
    opt = runtime= vector(mode="list",length=N)
  
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) = names(F1) = 
    names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA,nrep)
  
  # cat("checkpoint2:",p)
  
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xval, yval
    xy.obj = sim.xy(n,p,nval,rho,s,beta.type,snr)
    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        
        # Apply the regression method in hand
        if (indicator[j]) {
          runtime[[j]][i] = system.time({
            reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y,xy.obj$xval,xy.obj$yval)
          })[1]
        }else{
          runtime[[j]][i] = system.time({
            reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y)
          })[1]
        }
        
        
        # cat("checkpoint4.4 lasso out:",p)
        # cat("checkpoint4.5 lasso out:",coef(reg.obj))
        
        # Grab the estimated coefficients, and the predicted values on the
        # training and validation sets
        betahat = as.matrix(coef(reg.obj))
        m = ncol(betahat); nc = nrow(betahat)
        
        # cat("checkpoint5:",p)
        
        # Check for intercept
        if (nc == p+1) {
          intercept = TRUE
          betahat0 = betahat[1,]
          betahat = betahat[-1,]
        }else intercept = FALSE
        
        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =
            risk[[j]] = nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }
        
        
        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        delta = betahat - xy.obj$beta
        risk[[j]][i,] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i,] = risk[[j]][i,] + betahat0^2
        err.test[[j]][i,] = risk[[j]][i,] + xy.obj$sigma^2
        prop[[j]][i,] = 1 - err.test[[j]][i,] / err.null[i]
        nzs[[j]][i,] = colSums(betahat!=0)
        tpos = colSums((betahat!=0)*(xy.obj$beta!=0))
        fpos[[j]][i,] = nzs[[j]][i,]-tpos
        fneg[[j]][i,] = colSums((betahat==0)*(xy.obj$beta!=0))
        F1[[j]][i,] = 2*tpos/(2*tpos+fpos[[j]][i,]+fneg[[j]][i,])
        opt[[j]][i,] = (err.test[[j]][i,] - err.train[[j]][i,]) /
          err.train[[j]][i,]
        
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }
    
    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,fneg,F1,opt,sigma,runtime),file=file)
    }
  }
  
  
  # Save results now (in case of an error that might occur below)
  out = enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,nzs,fpos,
               fneg,F1,opt,sigma,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  out = choose.tuning.params(out)
  
  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}








