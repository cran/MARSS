## Wrapper function using new classes marssm, popWrap, marssMLE
## Default values are in MARSSsettings

MARSS = function(y,
    inits=NULL,
    constraint=NULL,
    fixed=NULL, free=NULL, 
    miss.value=-99,
    method = "kem",
    fit=TRUE, 
    silent = FALSE,
    control = NULL 
    ) 
{

  ## The popWrap() call does the following:  
  ## Set up default if some params left off
  ## Check that the user didn't pass in any illegal arguments (via call to checkMARSS())
  ## Set the initial conditions and make them the correct dimensions
  ## Translate Z constraint to fixed/free
  ## wrapperObj is an object of class "popWrap":
  ## list(data, m, inits, constraint, fixed, free, miss.value, control)

  if(!(method %in% allowed.methods)) stop(paste(method, "is not among the allowed methods. See ?MARSS."))
  
  this.method.allows = allowed[[method]]
  
  wrapperObj <- popWrap(y=y, this.method.allows, inits=inits, constraint=constraint, fixed=fixed, free=free, miss.value=miss.value, method=method, control=control, silent=silent)

  ## The as.marssm() call does the following: 
  ## Translate constraints names (shortcuts) to full fixed and free matrices
  ## modelObj is an obj of class "marssm":
  ## list(fixed, free, data, M, miss.value)

  modelObj = as.marssm(wrapperObj)
  
  ## Check that the model is ok
  tmp = is.marssm(modelObj)
    if(!isTRUE(tmp)) {
      if( !silent || silent==2 ) cat(tmp,"\n")
      stop("MARSS: marss model object is incomplete or inconsistent.\n", call.=FALSE)
    }

    X.names=NA
    #popWrap may have changed/set constraint$Z so need to use the wrapperObj$constraint$Z
    if(is.factor(wrapperObj$constraint$Z) ) X.names=unique(wrapperObj$constraint$Z)
    if( (identical(wrapperObj$constraint$Z,"use fixed/free") || is.matrix(wrapperObj$constraint$Z)) && !is.null(colnames(modelObj$fixed$Z)))
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
      if( !silent || silent==2 ) cat(tmp,"\n")
      cat("Stopped in MARSS: marssMLE object is incomplete or inconsistent.\n")
      cat("The incomplete/inconsistent MLE object is being returned. \n")
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
        }else if(MLEobj$convergence%in%c(10) && method %in% kem.methods){
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
