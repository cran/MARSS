#######################################################################################################
#   MARSSboot
#   This creates bootstrap parameter estimates
#   This is an MLE function and uses marssMLE objects
#   return(list(boot.params=boot.params, boot.data=boot.data, model=MLEobj, nboot=nboot, output=output, sim=sim, param.gen=param.gen, control=control))
#######################################################################################################
MARSSboot = function(MLEobj, nboot=1000, output="parameters", sim="parametric", 
  param.gen="KalmanEM", control=NULL, silent=FALSE) {

  # all equation numbers refer to Chapter 6, Shumway & Stoffer
  # tSteps is the number of time steps to do in each bootstrap of the data
  #     the default is to take this from the number of years in the original data
  # nboot is the number of bootstrap replicates to do
  # MLEobj is a marssMLE object (needs the structure and parameter elements of the list)
  # output tells what output to produce; this can be a vector
  #      data : returns the simulated (or bootstrapped) data.
  #      parameters : return the ML parameter estimates from the bootstrapped data
  #      all : return both
  # sim = the type of bootstrapping to do to create boot.data.  Choices are
  #      parametric uses parametric bootstrapping
  #      innovations uses Stoffer and Walls method
  # param.gen = how to generate the parameters.
  #      KalmanEM
  #      hessian
  # control is a list which holds options for the estimation function (see help file)
  # silent controls whether a progress bar is shown
  
  ###### Error-checking on the arguments
  if(!is.numeric(nboot))
    stop("MARSSboot: incorrect function arg. nboot must be numeric")
  if(any(nboot != trunc(nboot), nboot < 0, length(nboot)!=1) ) 
     stop("MARSSboot: incorrect function arg. nboot must be a positive non-zero integer")
  if(FALSE %in% (output %in% c("all","data","parameters"))) 
      stop("MARSSboot: incorrect function arg. output must be either all, data, parameters (in quotes)")
  if(FALSE %in% (sim %in% c("parametric","innovations"))) 
      stop("MARSSboot: incorrect function arg. sim must be either parametric or innovations (see help file)")
  if(FALSE %in% (param.gen %in% c("KalmanEM", "hessian"))) 
      stop("MARSSboot: incorrect function arg. est.method must be either KalmanEM or hessian (see help file)")
  if(FALSE %in% (silent %in% c(TRUE, FALSE))) 
      stop("MARSSboot: incorrect function arg. silent must be TRUE or FALSE")

  #mle.object$control is control args for estimation routine
  #control here is just if you want to change the control values in mle.object$control
  #code replaces mle.object$control element with the one in the passed in control obj
  if (!is.null(control)) {
    if(!is.list(control)) stop("MARSSboot: incorrect MARSSbootstrap arg. control must be a list.")
    EM.null = NULL
    en = c("minit","maxit","abstol","iter.V0")
    for (el in en) {
      null.flag = ( is.null(control[[el]]) )
 
      if (!null.flag) {
        null.flag = (!is.numeric(control[[el]]) || control[[el]] <= 0)	   
        if(!null.flag) MLEobj$control[[el]] = control[[el]]	#if control$el passed in and it is ok, reset MLEobj$control 
      }      
      EM.null <- c(EM.null, null.flag)  
    }  
    problem <- any(EM.null)
    if (problem) {
      if (!silent) cat("MARSSboot: Invalid control", en[EM.null], "\n")
      stop("MARSSboot: incorrect MARSSbootstrap arg -> control")
    }
  ##  trace must be set to (0) FALSE otherwise too much memory would be consumed
    MLEobj$control$trace = 0
    MLEobj$control$silent = TRUE  #prevent output from the estimation function
  }

  # Check for marssMLE properness; check that MLEobj$model$data exists and MLEobj$model$miss.value exists
  tmp = is.marssMLE(MLEobj)
  if(!isTRUE(tmp)) {
    if(!silent) print(tmp)
    stop("MARSSboot: MLE object is incomplete or inconsistent.")
  }
  # Check that it has par added on
  if(is.null(MLEobj$par)) stop("MARSSboot: MLE object is missing parameter estimates") 

  if(output=="all") output=c("data","parameters")
  
  # Check that if param.gen=hessian, the user hasn't set output to just "data"
  if(all(param.gen=="hessian", "data" %in% output))  {
    if(!MLEobj$control$silent) warning("MARSSboot: when param.gen=hessian, no boot.data are output")
    output = output[output!="data"]
    if(length(output)==0) stop("MARSSboot: when param.gen=hessian, output arg must be set to parameters")
    }
  if(param.gen=="hessian") sim = "none"

  #A little renamimg for code readability
  model=MLEobj$model 
  ##### Now check for any inconsistency in the passed in arguments
  if( (model$miss.value %in% model$data) & sim=="innovations") 
    stop("MARSSboot: innovations bootstrapping uses the innovations resampling and can only be done if there are no missing values in the data")

  ##### Set holders for output
  boot.params=NA; boot.data=NA  #dummy values 
  # these are the arrays for output
  if("data" %in% output) boot.data = array(NA,dim=c(dim(as.matrix(model$data)),nboot))
  if("parameters" %in% output) boot.params=array(NA,dim=c(length(MARSSvectorizeparam(MLEobj)),nboot))  
 
  ##### Read in model parameters
  m = dim(MLEobj$par$Z)[2]
  TT = dim(model$data)[2]  # length of time series
  n = dim(MLEobj$par$Z)[1]   

  ##### If using hessian to generate boot parameters, check if parSigma is already set
  if(param.gen=="hessian" & (is.null(MLEobj$parSigma) | is.null(MLEobj$parMean))){
    if(!silent) cat("MARSSboot: Computing the Hessian.  This might take awhile.\n")
    MLEobj = MARSShessian(MLEobj)  #adds the hessian; parSigma; and parMean matrix onto the mle.model object
    if(is.null(MLEobj$parSigma))  #if solve(hessian) didn't work in emHessian() then it sets parSigma to NULL
      stop("MARSSboot: Hessian could not be inverted to estimate the parameter var-cov matrix")
  }
		 
  ##### Set up the progress bar
  drawProgressBar = FALSE #If the time library is not installed, no prog bar
  if(!silent) { #then we can draw a progress bar
    prev = progressBar()
    drawProgressBar = TRUE
    }

  #####This part creates and stores the bootstrap data (not parameters yet just data)
  if(sim == "parametric")
    boot.data = MARSSsimulate(MLEobj$par, miss.loc=model$data, miss.value=model$miss.value, tSteps=TT, nsim=nboot )$sim.data  # make new data

  if(sim == "innovations")
    boot.data = MARSSinnovationsboot(MLEobj, nboot=nboot )$boot.data

  #####This part generates the bootstrap parameter
  if("parameters" %in% output){
    for(i in 1:nboot) {
      if( param.gen == "KalmanEM" ) {   
	        newmod = MLEobj$model
	        boot.control = MLEobj$control; boot.control$silent = TRUE #turn off output from MARSSkem()
	        newmod$data = array(boot.data[,,i], dim=dim(boot.data)[1:2])
          mle.object = list(model=newmod, start=MLEobj$start, control=boot.control )
	        boot.model = MARSSkem(mle.object)
	        boot.params[,i] = MARSSvectorizeparam(boot.model)
	# b495
	if (i == 1) rownames(boot.params) = names(MARSSvectorizeparam(boot.model))
      } #if KalmanEM

      if( param.gen == "hessian" ) {
	# b495
	if (i == 1) rownames(boot.params) = names(MLEobj$parMean)
        hess.params = rmvnorm(1, mean=MLEobj$parMean, sigma=MLEobj$parSigma, method="chol")
        boot.params[,i] = hess.params
      }   #if hessian   
    
      # Draw the progress bar if silent=F and time library is installed
    	if(drawProgressBar) prev <- progressBar(i/nboot,prev)
      }    #end nboot loop
    } # end if parameters in output
      
  return(list(boot.params=boot.params, boot.data=boot.data, model=MLEobj, nboot=nboot, output=output, sim=sim, param.gen=param.gen))
}
