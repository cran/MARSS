MARSS = function(y,
                 inits=NULL,
                 model=NULL,
                 miss.value=as.numeric(NA),
                 method = "kem",
                 form = "marxss",
                 fit=TRUE, 
                 silent = FALSE,
                 control = NULL,
                 MCbounds = NULL,
                 fun.kf = "MARSSkfas",
                 ... 
) 
{
  ## Start by checking the data, since if the data have major problems then the rest of the code
  ## will have problems
  if(!missing(miss.value)){
    stop("miss.value is deprecated in MARSS.  Replace missing values in y with NA.\n")
  }
  if(is.null(y)){ stop("MARSS: No data (y) passed in.",call.=FALSE) }
  if(is.ts(y)){ 
    stop("MARSS: Please convert the ts object to a matrix with time going across columns.\nThis command will work to convert the data: t(as.data.frame.ts(y)).\nThis command will make a row of the frequency info (if you need this as a covariate): t(as.data.frame.ts(stats:::cycle.ts(y))).\nType MARSSinfo(\"ts\") for more info.",call.=FALSE)
  }
  if(!(is.vector(y) | is.matrix(y))) stop("MARSS: Data (y) must be a vector or matrix (time going across columns).",call.=FALSE)
  if(length(y)==0) stop("MARSS: Data (y) is length 0.",call.=FALSE)
  if(is.vector(y)) y=matrix(y,nrow=1)
  if(any(is.nan(y))) cat("MARSS: NaNs in data are being replaced with NAs.  There might be a problem is NaNs shouldn't be in the data.\nNA is the normal missing value designation.\n")
  y[is.na(y)]=as.numeric(NA)
  
  MARSS.call = list(data=y, inits=inits, MCbounds=MCbounds, model=model, control=control, method=method, form=form, silent=silent, fit=fit, fun.kf=fun.kf, ...)
  
  #First make sure specified equation form has a corresponding function to do the conversion to marssMODEL (form=marss) object
  as.marss.fun = paste("MARSS.",form[1],sep="")
  tmp=try(exists(as.marss.fun,mode="function"),silent=TRUE)
  if(!isTRUE(tmp)){
    msg=paste(" MARSS.", form[1], "() function to construct a marssMODEL (form=marss) object does not exist.\n",sep="")
    cat("\n","Errors were caught in MARSS \n", msg, sep="") 
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call.=FALSE)
  }
  
  #Build the marssMODEL object from the model argument to MARSS()
  ## The as.marss.fun() call adds the following to MARSS.inputs list:
  ## $marss Translate model strucuture names (shortcuts) into a marssMODEL (form=marss) object put in $marss
  ## $alt.forms with any alternate marssMODEL objects in other forms that might be needed later
  ## a marssMODEL object is a list(data, fixed, free, tinitx, diffuse) 
  ## with attributes model.dims, X.names, form, equation
  ## error checking within the function is a good idea though not required
  ## if changes to the control values are wanted these can be set by changing MARSS.inputs$control
  if(silent==2) cat("Building the marssMODEL object from the model argument to MARSS().\n")
  MARSS.inputs = eval(call(as.marss.fun, MARSS.call))
  marss.object=MARSS.inputs$marss
  if(silent==2 ) cat(paste("Resulting model has ",attr(marss.object,"model.dims")$x[1]," state processes and ",attr(marss.object,"model.dims")$y[1]," observation processes.\n",sep=""))
  
  ## Check that the marssMODEL object output by MARSS.form() is ok
  ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
  if(silent==2) cat(paste("Checking that the marssMODEL object output by MARSS.",form,"() is ok.\n",sep=""))
  tmp = is.marssMODEL(marss.object)
  if(!isTRUE(tmp)) {
    if( !silent || silent==2 ) { cat(tmp); return(marss.object) } 
    stop("Stopped in MARSS() due to problem(s) with model specification. If silent=FALSE, marssMODEL object (in marss form) will be returned.\n", call.=FALSE)
  }
  
  ## checkMARSSInputs() call does the following:  
  ## Check that the user didn't pass in any illegal arguments
  ## and fill in defaults if some params left off
  ## This does not check model since the marssMODEL object is constructed
  ## via the MARSS.form() function above
  if(silent==2) cat("Running checkMARSSInputs().\n")
  MARSS.inputs=checkMARSSInputs(MARSS.inputs, silent=FALSE)
  
  
  ###########################################################################################################
  ##  MODEL FITTING
  ###########################################################################################################
  
  ## MLE estimation
  if(method %in% c(kem.methods, optim.methods)) {
    
    ## Create the marssMLE object
    
    MLEobj = list(marss=marss.object, model=MARSS.inputs$model, control=c(MARSS.inputs$control, list(MCbounds=MARSS.inputs$MCbounds), silent=silent), method=method, fun.kf=fun.kf)
    #Set the call form since that info needed for MARSSinits
    if(MLEobj$control$trace != -1){ MLEobj$call=MARSS.call } 
    
    # This is a helper function to set simple inits for a marss MLE model object
    MLEobj$start = MARSSinits(MLEobj, MARSS.inputs$inits)
    
    class(MLEobj) = "marssMLE"
    
    ## Check the marssMLE object
    ## is.marssMLE() calls is.marssMODEL() to check the model,
    ## then checks dimensions of initial value matrices.
    ## it also checks the control list and add defaults if some values are NULL
    tmp = is.marssMLE(MLEobj)
    
    #if errors, tmp will not be true, it will be error messages
    if(!isTRUE(tmp)) {
      if( !silent || silent==2) {
        cat(tmp)
        cat(" The incomplete/inconsistent MLE object is being returned.\n")
      }
      cat("Error: Stopped in MARSS() due to marssMLE object incomplete or inconsistent. \nPass in silent=FALSE to see the errors.\n\n")
      MLEobj$convergence=2
      return(MLEobj)
    }
    
    if(MLEobj$control$trace != -1){
      MLEobj.test=MLEobj
      MLEobj.test$par=MLEobj$start
      #Do this test with MARSSkfss since it has internal tests for problems
      kftest=try(MARSSkfss( MLEobj.test ), silent=TRUE)
      if(inherits(kftest, "try-error")){ 
        cat("Error: Stopped in MARSS() before fitting because MARSSkfss stopped.  Something is wrong with \n the model structure that prevents Kalman filter running.\n Try using control$trace=-1 and fit=FALSE and then look at the model that MARSS is trying to fit.\n")      
        MLEobj$convergence=2
        return(MLEobj.test)
      }
      if(!kftest$ok){
        cat(kftest$msg) 
        cat(paste("Error: Stopped in MARSS() before fitting because the MARSSkfss (Kalman filter) function returned errors.  Something is wrong with \n the model structure. Try using control$trace=-1 and fit=FALSE and then \n look at the model that MARSS is trying to fit.\n",kftest$errors,"\n",sep=""))
        MLEobj$convergence=2
        return(MLEobj.test)
      }
      MLEobj.test$kf=kftest
    }
    #Ey is needed for method=kem
    if(MLEobj$control$trace != -1 & MLEobj$method=="kem"){
      Eytest=try(MARSShatyt( MLEobj.test ), silent=TRUE)
      if(inherits(Eytest, "try-error")){ 
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt stopped.  Something is wrong with \n model structure that prevents MARSShatyt running.\n\n")
        MLEobj.test$convergence=2
        MLEobj.test$Ey=Eytest
        return(MLEobj.test)
      }
      if(!Eytest$ok){
        cat(Eytest$msg) 
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt returned errors.  Something is wrong with \n model structure that prevents function running.\n\n")      
        MLEobj.test$convergence=2
        MLEobj.test$Ey=Eytest
        return(MLEobj.test)
      }
      #MLEobj$Ey=Eytest
    }
    
    ## If a MCinit on the EM algorithm was requested
    if(MLEobj$control$MCInit) {
      MLEobj$start = MARSSmcinit(MLEobj)
    }
    
    if(!fit) MLEobj$convergence=3
    
    if(fit) {
      ## If no parameters are estimated, then set par element and get the states
      if(all(unlist(lapply(MLEobj$marss$free,is.fixed)))){
        MLEobj$convergence=3
        MLEobj$par=list(); for(el in attr(MLEobj$marss,"par.names")) MLEobj$par[[el]]=matrix(0,0,1)
        kf=MARSSkfss(MLEobj) #kfss is better here
        MLEobj$states=kf$xtT
        MLEobj$logLik=kf$logLik
        if(!is.null(kf[["VtT"]])){
          m = attr(MLEobj$marss,"model.dims")[["x"]][1]
          TT = attr(MLEobj$marss,"model.dims")[["data"]][2]
          if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
          if(m > 1) {
            states.se = matrix(0, nrow=m, ncol=TT)
            for(i in 1:TT) 
              states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
          }
        }else{  states.se=NULL }
        MLEobj[["states.se"]] = states.se
        Ey=MARSShatyt(MLEobj)
        MLEobj$ytT=Ey[["ytT"]]
        if(!is.null(Ey[["OtT"]])){
          n = attr(MLEobj$marss,"model.dims")[["y"]][1]
          TT = attr(MLEobj$marss,"model.dims")[["data"]][2]
          if(n == 1) y.se = sqrt(matrix(Ey[["OtT"]][,,1:TT], nrow=1))
          if(n > 1) {
            y.se = matrix(0, nrow=n, ncol=TT)
            for(i in 1:TT) 
              y.se[,i] = t(sqrt(takediag(Ey[["OtT"]][,,i])))
          }
        }else{  y.se=NULL }
        MLEobj$y.se=y.se
        if(MLEobj$control$trace>0){ 
          MLEobj$kf=MARSSkf(MLEobj)
          MLEobj$Ey=MARSShatyt(MLEobj)
        }
      }else{ #there is something to estimate
        if(silent==2 ) cat(paste("Fitting model with ",method,".\n",sep=""))
        ## Fit and add param estimates to the object
        if(method %in% kem.methods) MLEobj = MARSSkem(MLEobj)
        if(method %in% optim.methods) MLEobj = MARSSoptim(MLEobj)
      }
      
      #apply X and Y names various X and Y related elements
      MLEobj = MARSSapplynames(MLEobj)
      
      ## Add AIC and AICc and coef to the object
      ## Return as long as there are no errors, but might not be converged
      if((MLEobj$convergence%in%c(0,1)) | (MLEobj$convergence%in%c(10,11) && method %in% kem.methods) ){
        MLEobj = MARSSaic(MLEobj)
        MLEobj$coef = coef(MLEobj,type="vector")
        #add on the states.se's
        kf=MARSSkfss(MLEobj) #kfss is better here
        if(!is.null(kf[["VtT"]])){
          m = attr(MLEobj$marss,"model.dims")[["x"]][1]
          TT = attr(MLEobj$marss,"model.dims")[["data"]][2]
          if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
          if(m > 1) {
            states.se = matrix(0, nrow=m, ncol=TT)
            for(i in 1:TT) 
              states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
          }
        }else{  states.se=NULL }
         MLEobj[["states.se"]] = states.se
        }
    } # fit the model
    
    if((!silent || silent==2) & MLEobj$convergence %in% c(0,1,3,10,11,12)){ print(MLEobj) }
    if((!silent || silent==2) & !(MLEobj$convergence %in% c(0,1,3,10,11,12))){ cat(MLEobj$errors) } # 3 added since don't print if fit=FALSE
    if((!silent || silent==2) & !fit) print(MLEobj$model)
    
    return(MLEobj)
  } # end MLE methods
  
  return("method allowed but it's not in kem.methods or optim.methods so marssMLE object was not created")
}  

