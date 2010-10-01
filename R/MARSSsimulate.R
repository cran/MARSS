#######################################################################################################
#   MARSSsimulate function
#   Parametrically simulates from a MARSS parameter list
#######################################################################################################
MARSSsimulate = function(parList, tSteps=100, nsim=1, silent=TRUE, miss.loc=NULL, miss.value=NULL) {
  #parList is a list of the parameters as would come out of passing in MLEmodel$par
  # tSteps is the number of time steps to do in each bootstrap of the data
  # miss.loc is an optional (n x tSteps x nsim) matrix specifying where to put missing values
  # if miss.loc is the same for all nsim, can pass in dim=c(n, tSteps)
  if(is.null(parList$Z) || is.null(dim(parList$Z)))     
    stop("Stopped in MARSSsimulate() due to problem with Z matrix of the parameter list. Missing or it is not a matrix.\n", call.=FALSE)  
  n = dim(parList$Z)[1]
  m = dim(parList$Z)[2]
  
  ###### Error-checking on parameters
  tmp = MARSScheckpar(parList, n, m)
  if(!isTRUE(tmp)) {
    if(!silent) cat(tmp)
    stop("Stopped in MARSSsimulate() due to problem with parameter matrices.\n", call.=FALSE)
  }

  ###### Error-checking on the arguments
  msg=NULL
  if(!is.numeric(tSteps) || !is.numeric(nsim))
    msg=c(msg,"Non-numeric tSteps or nsim argument(s).\n")
  if(is.numeric(tSteps) && (tSteps != trunc(tSteps) || tSteps <= 0)) 
    msg=c(msg,"tSteps must be a positive non-zero integer.\n")
  if(is.numeric(nsim) && (nsim != trunc(nsim) || nsim <= 0))
    msg=c(msg,"nsim must be a positive non-zero integer.\n")
  if(!is.null(miss.loc) && (!isTRUE( all.equal(dim(miss.loc), c(n,tSteps)) ) 
    && !isTRUE(all.equal( dim(miss.loc), c(n,tSteps,nsim)) ) )  )
    msg=c(msg,"Incorrect input arg: miss.loc dim must be n x tSteps or n x tSteps x nsim, where n=dim(Z)[1].\n")
  if(!is.null(miss.loc) && is.null(miss.value)) 
    msg=c(msg,"If miss.loc is passed in, you must specify miss.value.\n")
  if(!is.null(miss.value) && !is.na(miss.value) && !is.numeric(miss.value)) 
    msg=c(msg,"miss.value must be numeric (or NA).\n")
  if(!is.null(miss.value) && length(miss.value)!=1) 
    msg=c(msg,"miss.value must be length 1.\n")
  if(!is.null(msg)){
    cat("\nErrors were caught in MARSSsimulate\n", msg) 
    stop("Stopped in MARSSsimulate() due to argument problem(s).\n", call.=FALSE)
  }
  if( !is.null(miss.loc) && identical(NA,miss.value) ){ miss.loc.TF = (is.na(miss.loc) & !is.nan(miss.loc)) }
  if( !is.null(miss.loc) && identical(NaN,miss.value) ){ miss.loc.TF = is.nan(miss.loc) }
  if( !is.null(miss.loc) && !is.na(miss.value) ){ miss.loc.TF = (miss.loc==miss.value) }
       
  ##### Set holders for output
  #if user passed in miss.loc dim=c(n,tSteps) assume they wanted that repeated for all nsim's
  #set up miss.loc if user didn't pass it in; 
  if(is.null(miss.loc))     #means no missing values
    miss.loc.TF=array(TRUE,dim=c(n,tSteps))  #1 means not missing
  if(length(dim(miss.loc.TF))==2) miss.loc.TF = array(miss.loc.TF, dim=c(n,tSteps,nsim))
  #sim.data = array(NA,dim=c(tSteps,n,nsim))
  #sim.states = array(NA,dim=c(tSteps,m,nsim))
  sim.data = array(NA,dim=c(n, tSteps, nsim))
  sim.states = array(NA,dim=c(m, tSteps, nsim))
 		 
  ##### Set up the progress bar
  drawProgressBar = FALSE #If the time library is not installed, no prog bar
  if(!silent) { #then we can draw a progress bar
    prev = progressBar()
    drawProgressBar = TRUE
    }

  ##### Set up holders
  newData = matrix(NA, n,tSteps)
  newStates = matrix(NA, m, tSteps+1) # States = years x subpops
  
  for( i in 1:nsim){
    newStates[,1] = parList$x0 # t = 0
    # create a matrices for observation error
    obs.error = t(rmvnorm(tSteps, mean = rep(0, n), sigma = parList$R, method="chol"))
    pro.error = t(rmvnorm(tSteps, mean = rep(0, m), sigma = parList$Q, method="chol"))
    for(j in 2:(tSteps+1)) {
      newStates[,j] = parList$B %*% newStates[,j-1] + parList$U + pro.error[,j-1]
      newData[,j-1] = parList$Z %*% newStates[,j] + parList$A + obs.error[,j-1]
    }
    if(!is.null(miss.value)) newData[miss.loc.TF[,,i]] = miss.value
    newStates=newStates[,2:(tSteps+1)] 
    sim.data[,,i] = as.matrix(newData)
    sim.states[,,i] = as.matrix(newStates)
    # reset newStates to its original dim
    newStates = matrix(NA, m, tSteps+1)
    # Draw the progress bar if silent=F and time library is installed
    if(drawProgressBar) prev <- progressBar(i/nsim,prev)
  } # end of for loop for nsim 
       
  return(list(sim.states=sim.states, sim.data=sim.data, par=parList, miss.loc=miss.loc, miss.value=miss.value, tSteps=tSteps, nsim=nsim))
}
