#######################################################################################################
#   Parameter estimation using R's optim function
#   Minimal error checking is done.  You should run is.marssMLE() before calling this.
#   Likelihood computation is via  Kalman filter
#######################################################################################################
MARSSoptim = function(MLEobj) {
# This function does not check if user specified a legal MLE object.
# free and fixed are list of constraint matrices.  Values that are not fixed must be designated NA in "fixed"; 
# Values that are not free must be designated NA in "free"
  model.el = names(MLEobj$model$free)[names(MLEobj$model$free)!="V0"]
  model.el.w.V0 = names(MLEobj$model$free)
  tmp = is.marssMLE(MLEobj)
  if(!isTRUE(tmp)) {
      cat(tmp)
      stop("Stopped in MARSSoptim() because marssMLE object is incomplete or inconsistent.\n", call.=FALSE)
    }
  tmp = describe.marssm(MLEobj$model)  
  for(elem in c("Q","R")){
    varcov.type = substr(tmp[[elem]],1,6)
    if( !is.fixed(MLEobj$model$fixed[[elem]]) && varcov.type!="diagon" && varcov.type!="scalar") {
       stop(paste("Stopped in MARSSoptim(). If estimated, ",elem," must be diagonal.\n", sep=""), call.=FALSE)  }
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
  
  ## Disallow V0 estimation
  if(!is.fixed(fixed$V0)) stop("Stopped in MARSSoptim(). V0 cannot be estimated.\n", call.=FALSE)
       
  # If V0 is fixed to be zero, then the Newton algorithms are run with a diag V0 set large.  At end, the kalman filter is rerun with 
  if(!is.fixed(fixed$x0)){
   if(!identical(unname(fixed$V0), array(0,dim=c(m,m)))){ 
      stop("Stopped in MARSSoptim(). If x0 is estimated, V0 must be 0.  See discussion regarding initial conditions in manual.\n",call.=FALSE)
   }else{
        D=as.design(fixed$x0, free$x0)$D  #need this many places
        V0 = control$iter.V0 * D%*%t(D) #if some x0 are shared, they need V0 with 100% correlation
        }
   }else { V0 = fixed$V0 }   #if x0 is fixed, x0 is treated as a prior
  tmp.inits$V0=V0 #set to whatever V0 is fixed to; again for returning; not used

    
  #The code is used to set things up to use MARSSvectorizeparam to just select inits for the estimated parameters
  tmp.MLEobj = MLEobj
  tmp.MLEobj$model$fixed$V0 = V0 #set up the V0 used in optim search
  tmp.MLEobj$par = tmp.inits  #set initial conditions for estimated parameters
  #diagonal elements of Q and R will be logged' need to log both par and fixed
  diag(tmp.MLEobj$par$Q) = log(diag(tmp.MLEobj$par$Q))
  diag(tmp.MLEobj$par$R) = log(diag(tmp.MLEobj$par$R))
  diag(tmp.MLEobj$model$fixed$Q) = log(diag(tmp.MLEobj$model$fixed$Q))
  diag(tmp.MLEobj$model$fixed$R) = log(diag(tmp.MLEobj$model$fixed$R))
  # will return the inits only for the estimated parameters and diag elements of Q and R will be logged
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
      diag(tmp.MLEobj$par$Q) = exp(diag(tmp.MLEobj$par$Q))
      diag(tmp.MLEobj$par$R) = exp(diag(tmp.MLEobj$par$R))
      pars = MARSSvectorizeparam(tmp.MLEobj)
      #now put the estimated values back into MLEobj; fixed values set by $fixed
      #so V0 is back to zero if it had been set to a temporary value for estimation
      MLEobj.return = MARSSvectorizeparam(MLEobj.return, pars)
      kf = MARSSkf(MLEobj.return$model$data, MLEobj.return$par, miss.value = MLEobj.return$model$miss.value, init.state="x10")
      }else{
      if(optim.output$convergence==10) optim.output$message=c("degeneracy of the Nelder-Mead simplex\n",optim.output$message)
      optim.output$counts = NULL      
      if( !control$silent ) cat("MARSSoptim() stopped with errors. No parameter estimates returned.\n")
      if( control$silent==2 ) cat("MARSSoptim() stopped with errors. No parameter estimates returned. See $errors in output for details.\n")
 
      MLEobj.return$par = NULL
      MLEobj.return$errors = optim.output$message
      kf = NULL
      }
  
  if(!is.null(kf)){
    MLEobj.return$kf = kf
    MLEobj.return$states = kf$xtT
    MLEobj.return$numIter = optim.output$counts[1]
    MLEobj.return$logLik = kf$logLik
  }
  MLEobj.return$method = MLEobj$method
  
  ## Add AIC and AICc to the object
  if(!is.null(kf)) MLEobj.return = MARSSaic(MLEobj.return)

  ## Calculate confidence intervals based on state std errors, see caption of Fig 6.3 (p337) Shumway and Stoffer
  if(!is.null(kf)){
    TT = dim(MLEobj.return$model$data)[2]; m = dim(MLEobj.return$model$fixed$Q)[1]
    if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
    if(m > 1) {
      states.se = matrix(0, nrow=m, ncol=TT)
      for(i in 1:TT) states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
    }
    MLEobj.return$states.se = states.se
    }
    
  return(MLEobj.return)
}

neglogLik = function(x, MLEobj=NULL){  #NULL assignment needed for optim call syntax
    MLEobj = MARSSvectorizeparam(MLEobj, x)
    diag(MLEobj$par$Q) = exp(diag(MLEobj$par$Q))
    diag(MLEobj$par$R) = exp(diag(MLEobj$par$R))
    negLL = MARSSkf(MLEobj$model$data, MLEobj$par, miss.value = MLEobj$model$miss.value, init.state="x10")$logLik
    -1*negLL
     }
