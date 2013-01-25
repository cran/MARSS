MARSS = function(y,
    inits=NULL,
    model=NULL,
    miss.value=NA,
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
if(is.na(miss.value)) miss.value=as.numeric(miss.value)

MARSS.call = list(data=y, inits=inits, MCbounds=MCbounds, model=model, miss.value=miss.value, control=control, method=method, form=form, silent=silent, fit=fit, fun.kf=fun.kf, ...)

#First make sure specified equation form has a corresponding function to do the conversion to marssm object
  as.marssm.fun = paste("MARSS.",form,sep="")
  tmp=try(exists(as.marssm.fun,mode="function"),silent=TRUE)
  if(!isTRUE(tmp)){
    msg=paste(" MARSS.", form, "() function to convert form to marss object does not exist.\n",sep="")
    cat("\n","Errors were caught in MARSS \n", msg, sep="") 
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call.=FALSE)
  }

#Build the marssm object from the model argument to MARSS()
  ## The as.marssm.fun() call does the following: 
  ## Translate model strucuture names (shortcuts) into a marssm object
  ## which is added to the MARSS.inputs list
  ## is is a list(data, fixed, free, miss.value, X.names, tinitx, diffuse)
  ## error checking within the function is a good idea though not required
  ## if changes to the control values are wanted these can be set by changing MARSS.inputs$control
  ## Add list element form.info if there in information about the form (say for printing) that you want access to later
  as.marssm.fun = paste("MARSS.",form,sep="")
  MARSS.inputs = eval(call(as.marssm.fun, MARSS.call))
  modelObj=MARSS.inputs$marssm
  
  ## Check that the marssm object output by MARSS.form() is ok
  ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
  tmp = is.marssm(modelObj)
    if(!isTRUE(tmp)) {
      if( !silent || silent==2 ) { cat(tmp); return(modelObj) } 
      stop("Stopped in MARSS() due to problem(s) with model specification. If silent=FALSE, modelObj will be returned.\n", call.=FALSE)
    }

  #apply some generic row and col naming to fixed and free matrices in modelObj
  modelObj = MARSSapplynames(modelObj, X.names=modelObj$X.names)

  ## checkMARSSInputs() call does the following:  
  ## Check that the user didn't pass in any illegal arguments
  ## and fill in defaults if some params left off
  ## This does not check model since the marssm object is constructed
  ## via the MARSS.form() function above
  MARSS.inputs=checkMARSSInputs(MARSS.inputs, silent=FALSE)


###########################################################################################################
##  MODEL FITTING
###########################################################################################################

  ## MLE estimation
  if(method %in% c(kem.methods, optim.methods)) {
    
    ## Create the marssMLE object

    MLEobj = list(model=modelObj, control=c(MARSS.inputs$control, list(MCbounds=MARSS.inputs$MCbounds), silent=silent), method=method, fun.kf=fun.kf, form=form)
    MLEobj$form.info=MARSS.inputs$form.info
    # This is a helper function to set simple inits for a marss MLE model object
    MLEobj$start = MARSSinits(MLEobj, MARSS.inputs$inits)

    class(MLEobj) = "marssMLE"

    ## Check the marssMLE object
    ## is.marssMLE() calls is.marssm() to check the model,
    ## then checks dimensions of initial value matrices.
    ## it also checks the control list and add defaults if some values are NULL
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
      ## If not parameters are estimated, then get the states
      if(all(unlist(lapply(MLEobj$model$free,is.fixed)))){
        MLEobj$convergence=3
        MLEobj$par=list(); for(el in names(MLEobj$model$free)) MLEobj$par[[el]]=matrix(0,0,1)
        kf=MARSSkf(MLEobj) 
        MLEobj$states=kf$xtT
        MLEobj$logLik=kf$logLik
        if(!is.null(kf[["VtT"]])){
          m = dim(MLEobj$model$fixed$x0)[1]
          TT = dim(MLEobj$model$data)[2]
          if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
          if(m > 1) {
          states.se = matrix(0, nrow=m, ncol=TT)
          for(i in 1:TT) 
            states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
          }
        }else  states.se=NULL
        MLEobj$states.se = states.se
        MLEobj$ytT=MARSShatyt(MLEobj)$ytT
        if(MLEobj$control$trace>0){ 
          MLEobj$kf=MARSSkf(MLEobj)
          MLEobj$Ey=MARSShatyt(MLEobj)
        }
      }else{ #there is something to estimate
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
      }
      } # fit the model
  
   if(MLEobj$control$trace == -1){ # just save the form information from $call
     MLEobj$call$form = MARSS.call$form
   }else{ MLEobj$call=MARSS.call } 
      
   if(!silent & MLEobj$convergence %in% c(0,1,3,10,11)){ print(MLEobj) }
   if(!silent & !(MLEobj$convergence %in% c(0,1,3,10,11))){ cat(MLEobj$errors) } # 3 added since don't print if fit=FALSE
   if(!silent & !fit) print(modelObj)

   return(MLEobj)
  } # end MLE methods

  return("method allowed but it's not in kem.methods or optim.methods so marssMLE object was not created")
}  

