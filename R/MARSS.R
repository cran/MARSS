## Wrapper function using new classes marssm, popWrap, marssMLE
## Default values are in MARSSsettings

MARSS = function(y,
    inits=NULL,
    model=NULL,
    miss.value=NA,
    method = "kem",
    fit=TRUE, 
    silent = FALSE,
    control = NULL 
    ) 
{
  fixed=NULL; free=NULL

  ## The popWrap() call does the following:  
  ## Set up default if some params left off
  ## Check that the user didn't pass in any illegal arguments (via call to checkMARSS())
  ## Set the initial conditions and make them the correct dimensions
  ## Translate Z model structure to fixed/free
  ## wrapperObj is an object of class "popWrap":
  ## list(data, m, inits, model, fixed, free, miss.value, control)

  if(!(method %in% allowed.methods)){
    msg=paste(" ", method, "is not among the allowed methods. See ?MARSS.\n")
    cat("\n","Errors were caught in MARSS \n", msg, sep="") 
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call.=FALSE)
  }
  
  this.method.allows = allowed[[method]]
  
  wrapperObj <- popWrap(y=y, this.method.allows, inits=inits, model=model, fixed=fixed, free=free, miss.value=miss.value, method=method, control=control, silent=silent)

  ## The as.marssm() call does the following: 
  ## Translate model strucuture names (shortcuts) to full fixed and free matrices
  ## modelObj is an obj of class "marssm":
  ## list(fixed, free, data, M, miss.value)

  modelObj = as.marssm(wrapperObj)
  
  ## Check that the model is ok
  tmp = is.marssm(modelObj)
    if(!isTRUE(tmp)) {
      if( !silent || silent==2 ) cat(tmp) 
      stop("Stopped in MARSS() due to problem(s) with model specification.\n", call.=FALSE)
    }

    X.names=NA
    #popWrap may have changed/set model$Z so need to use the wrapperObj$model$Z
    if(is.factor(wrapperObj$model$Z) ) X.names=unique(wrapperObj$model$Z)
    if( is.matrix(wrapperObj$model$Z) && !is.null(colnames(modelObj$fixed$Z)))
      X.names=colnames(modelObj$fixed$Z)
    #apply some generic naming
    modelObj = MARSSapplynames(modelObj, X.names=X.names)

###########################################################################################################
##  MODEL FITTING
##  THE REST OF THIS FRONTEND IS NOT MARSS() SPECIFIC
###########################################################################################################

  ## MLE estimation
  if(method %in% c(kem.methods, optim.methods)) {
    
    ## Create the marssMLE object
    # This is a helper function to set simple inits for a marss MLE model object
    wrapperObj$inits = MARSSinits(modelObj, wrapperObj$inits)

    MLEobj = list(model=modelObj, start=wrapperObj$inits, control=c(wrapperObj$control, silent=silent), method=method)
    class(MLEobj) = "marssMLE"

    ## Check the marssMLE object
    ## is.marssMLE() calls is.marssm() to check the model,
    ## then checks dimensions of initial value matrices.
    tmp = is.marssMLE(MLEobj)
    if(!isTRUE(tmp)) {
      if( !silent ) {
        cat(tmp)
        cat(" The incomplete/inconsistent MLE object is being returned.\n")
        }
      cat("Error: Stopped in MARSS() due to marssMLE object incomplete or inconsistent. \nPass in silent=FALSE to see the errors.\n\n")
      MLEobj$convergence=2
      return(MLEobj)
    }

    ## If a MCinit on the EM algorithm was requested
    if(MLEobj$control$MCInit) {
      MLEobj$start = MARSSmcinit(MLEobj)
    }
    
    if(!silent & !fit) print(modelObj)

    if(fit) {
      ## Fit via EM and add param estimates to the object 
      if(method %in% kem.methods) MLEobj = MARSSkem(MLEobj)
      if(method %in% optim.methods) MLEobj = MARSSoptim(MLEobj)
    
      ## Add AIC and AICc to the object
      ## Return as long as there are no errors, but might not be converged
      if(MLEobj$convergence%in%c(0,1)) {
        MLEobj = MARSSaic(MLEobj)
        if(!silent){ print(MLEobj) }
        }else if(MLEobj$convergence%in%c(3,10,11) && method %in% kem.methods){
          MLEobj = MARSSaic(MLEobj)
          if(!silent ){ print(MLEobj) }
        }else if(!silent) cat(MLEobj$errors)  #stopped with errors
      }
      #apply names to the start and par elements
      MLEobj = MARSSapplynames(MLEobj)

      return(MLEobj)

  } # end MLE methods

  return("method allowed but it's not in kem.methods or optim.methods so marssMLE object was not created")
}  

