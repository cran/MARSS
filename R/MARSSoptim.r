#######################################################################################################
#   Parameter estimation using R's optim function
#   Minimal error checking is done.  You should run is.marssMLE() before calling this.
#   Likelihood computation is via  Kalman filter
#######################################################################################################
MARSSoptim = function(MLEobj) {
# This function does not check if user specified a legal MLE object.
# free and fixed are a list of model matrices.  Values that are not fixed must be designated NA in "fixed"; 
# Values that are not free must be designated NA in "free"
  tmp = is.marssMLE(MLEobj)
  if(!isTRUE(tmp)) {
      cat(tmp)
      stop("Stopped in MARSSoptim() because marssMLE object is incomplete or inconsistent.\n", call.=FALSE)
    }
  
  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MLEobj$model$data #must have time going across columns
  model = MLEobj$model
  free = MLEobj$model$free
  fixed = MLEobj$model$fixed
  tmp.inits = MLEobj$start
  control=MLEobj$control
  m=dim(free$Q)[1]; n=dim(free$R)[1]
  
  ## Set up the control list for optim; only pass in optim control elements
  control.names=c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr","pgtol", "temp", "tmax")
  optim.control=list()
  for(elem in control.names){
  if(!is.null(control[[elem]])) optim.control[[elem]]=control[[elem]]
  }
  if(is.null(control$lower)){ lower = -Inf
  }else lower=control$lower
  if(is.null(control$upper)){ upper = Inf
  }else upper=control$upper
    
  #The code is used to set things up to use MARSSvectorizeparam to just select inits for the estimated parameters
  tmp.MLEobj = MLEobj
  tmp.MLEobj$par = tmp.inits  #set initial conditions for estimated parameters
  #diagonal elements of Q and R will be logged' need to log both par and fixed
  for(elem in c("Q","R","V0")){
        the.par=tmp.MLEobj$par[[elem]]
        is.zero=diag(tmp.MLEobj$par[[elem]])==0
        if(any(is.zero)) diag(the.par)[is.zero]=1    #so the chol doesn't fail if there are zeros on the diagonal
        tmp.MLEobj$par[[elem]]=t(chol(the.par))
        if(any(is.zero)) diag(tmp.MLEobj$par[[elem]])[is.zero]=0
        tmp.MLEobj$model$free[[elem]][upper.tri(tmp.MLEobj$par[[elem]])]=NA
        is.na.par=is.na(diag(fixed[[elem]]))
        the.par= fixed[[elem]][!is.na.par,!is.na.par,drop=FALSE]
        if(!all(is.na.par)){
          is.zero= (diag(the.par)==0)
          if(any(is.zero)) diag(the.par)[is.zero]=1    #so the chol doesn't fail if there are zeros on the diagonal
          chol.the.par=t(chol(the.par))
          if(any(is.zero)) diag(chol.the.par)[is.zero]=0    
          tmp.MLEobj$model$fixed[[elem]][!is.na.par,!is.na.par]=chol.the.par
        }
    }
  # will return the inits only for the estimated parameters
  pars = MARSSvectorizeparam(tmp.MLEobj) 

  optim.output = try(optim(pars, neglogLik, MLEobj=tmp.MLEobj, method = MLEobj$method, lower = lower, upper = upper, control = optim.control, hessian = FALSE), silent=TRUE  )
  if(class(optim.output)=="try-error") optim.output = list(convergence=53, message=c(" MARSSkf() call used to compute log likelihood encountered numerical problems\n and could not return logLik. Sometimes better initial conditions helps.\n"))
  MLEobj.return=list(); class(MLEobj.return) = "marssMLE"
  MLEobj.return$iter.record=optim.output$message
  MLEobj.return$control=MLEobj$control
  MLEobj.return$model=MLEobj$model
  MLEobj.return$start = tmp.inits #set to what was used here
  MLEobj.return$convergence = optim.output$convergence
  if(optim.output$convergence %in% c(1,0)) {
      if((!control$silent || control$silent==2) && optim.output$convergence==0) cat(paste("Success! Converged in ",optim.output$counts[1]," interations.\n",sep=""))
      if((!control$silent || control$silent==2) && optim.output$convergence==1) cat(paste("Warning! Max iterations of ", control$maxit," reached before convergence.\n",sep=""))
      tmp.MLEobj = MARSSvectorizeparam(tmp.MLEobj, optim.output$par)
      #par has the fixed and estimated values with diags of Q and R logged
  for(elem in c("Q","R","V0")){
        L=tmp.MLEobj$par[[elem]]
        tmp.MLEobj$par[[elem]]=L%*%t(L)
    } #end for
    
      pars = MARSSvectorizeparam(tmp.MLEobj)  #now the pars values are adjusted back to normal scaling
      #now put the estimated values back into MLEobj; fixed values set by $fixed
      MLEobj.return = MARSSvectorizeparam(MLEobj.return, pars)
      kf.out = MARSSkf(MLEobj.return$model$data, MLEobj.return$par, miss.value = MLEobj.return$model$miss.value, init.state=control$kf.x0)
      }else{
      if(optim.output$convergence==10) optim.output$message=c("degeneracy of the Nelder-Mead simplex\n",optim.output$message)
      optim.output$counts = NULL      
      if( !control$silent ) cat("MARSSoptim() stopped with errors. No parameter estimates returned.\n")
      if( control$silent==2 ) cat("MARSSoptim() stopped with errors. No parameter estimates returned. See $errors in output for details.\n")
 
      MLEobj.return$par = NULL
      MLEobj.return$errors = optim.output$message
      kf.out = NULL
      }
  
  if(!is.null(kf.out)){
    MLEobj.return$kf = kf.out
    MLEobj.return$states = kf.out$xtT
    MLEobj.return$numIter = optim.output$counts[1]
    MLEobj.return$logLik = kf.out$logLik
  }
  MLEobj.return$method = MLEobj$method
  
  ## Add AIC and AICc to the object
  if(!is.null(kf.out)) MLEobj.return = MARSSaic(MLEobj.return)

  ## Calculate confidence intervals based on state std errors, see caption of Fig 6.3 (p337) Shumway and Stoffer
  if(!is.null(kf.out)){
    TT = dim(MLEobj.return$model$data)[2]; m = dim(MLEobj.return$model$fixed$Q)[1]
    if(m == 1) states.se = sqrt(matrix(kf.out$VtT[,,1:TT], nrow=1))
    if(m > 1) {
      states.se = matrix(0, nrow=m, ncol=TT)
      for(i in 1:TT) states.se[,i] = t(sqrt(takediag(kf.out$VtT[,,i])))
    }
    MLEobj.return$states.se = states.se
    }
    
  return(MLEobj.return)
}

neglogLik = function(x, MLEobj=NULL){  #NULL assignment needed for optim call syntax
    MLEobj = MARSSvectorizeparam(MLEobj, x)
  for(elem in c("Q","R","V0")){
        L=MLEobj$par[[elem]]
        MLEobj$par[[elem]]=L%*%t(L)
    } #end for
    if( MLEobj$control$kf.x0 == "x00"){
    negLL = MARSSkf(MLEobj$model$data, MLEobj$par, miss.value = MLEobj$model$miss.value, init.state=MLEobj$control$kf.x0 )$logLik    
    }else{ #must be x10
    negLL = MARSSkfas(MLEobj$model$data, MLEobj$par, miss.value = MLEobj$model$miss.value, init.state=MLEobj$control$kf.x0, diffuse=MLEobj$control$diffuse )$logLik
    }
    -1*negLL
     }
