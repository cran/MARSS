#######################################################################################################
#   MARSSmcinit function
#   Does a simple MonteCarlo initialization of the EM routine
#######################################################################################################

MARSSmcinit = function(MLEobj) {

  drawProgressBar = FALSE #If the time library is not installed, no prog bar
  if(!MLEobj$control$silent) { #then we can draw a progress bar
    cat("\n"); cat("> Starting Monte Carlo Initializations\n")
    prev = progressBar()
    drawProgressBar = TRUE
  }

  y = MLEobj$model$data
  ## M matrix for handling missing values
  M = MLEobj$model$M
  #Make sure the missing vals in y are zeroed out
  for(i in 1:dim(y)[2]){ y[!as.logical(takediag(M[,,i])),i]=0 }

  m = dim(MLEobj$model$fixed$Z)[2]
  n = dim(y)[1]
  TT = dim(y)[2]
  free.tmp = MLEobj$model$free 
  dim.tmp = list(Z=c(n,m), A=c(n,1), R=c(n,n), B=c(m,m), U=c(m,1), Q=c(m,m))
  # since boundsInits names are different
  tmp = MLEobj$control$boundsInits
  bounds.tmp = list(B=tmp$B, U=tmp$U, Q=tmp$logQ, Z=tmp$Z, A=tmp$A, R=tmp$logR)
  init = bestinits = MLEobj$start
  bestLL = -1.0e10

  # loop over numInits: # of random draws of initial values
  for(loop in 1:MLEobj$control$numInits) {
    init.loop = init
      
    # Draw random values
    en = c("Z", "A", "R", "B", "U", "Q")
    for(el in en) {
      if(FALSE %in% is.na(free.tmp[[el]])){
        dim.param = dim.tmp[[el]]
        bounds.param = bounds.tmp[[el]]
        tmp = table(free.tmp[[el]])
        free.levels = names(tmp)
        numGroups <- length(free.levels)
        Ztmp <- matrix(0, dim.param[1]*dim.param[2], numGroups)  # matrix to allow shared values
        for(i in free.levels) 
          Ztmp[which(as.vector(free.tmp[[el]])==i), which(free.levels==i)] <- 1
	if (el %in% c("Q", "R")) {
          element.random = array(exp(runif(numGroups, bounds.param[1], bounds.param[2])), dim=c(numGroups,1))
	}
	else {
	  element.random = array(runif(numGroups, bounds.param[1], bounds.param[2]), dim=c(numGroups,1))
	}
        param.random = array(Ztmp%*%element.random, dim = dim.param)
      }
      else param.random=0
   
      fix.tmp = MLEobj$model$fixed[[el]]
      fix.tmp[is.na(fix.tmp)] = 0 
      init.loop[[el]] = fix.tmp + param.random 
    }

    ## x0
    x0init = MLEobj$start$x0
    x.lo = ifelse(x0init > 0, 0.75*x0init, 1.25*x0init)
    x.hi = ifelse(x0init > 0, 1.25*x0init, 0.75*x0init)
    fix.tmp = MLEobj$model$fixed$x0
    fix.tmp[is.na(fix.tmp)] = 0 
    init.loop$x0 = fix.tmp + array(runif(m, x.lo, x.hi), dim=c(m,1))    

    ## Call MARSSkem() with these inits 
    MLEobj$start = init.loop
    MLEobj$control$maxit = MLEobj$control$numInitSteps
    MLEobj$control$minit = 1
    MLEobj$control$silent = TRUE #don't print convergence information during kem call          
    MLEobj = MARSSkem(MLEobj)    

    if(drawProgressBar==TRUE) prev = progressBar(loop/MLEobj$control$numInits, prev)

    ## Check whether the likelihood is the best observed
## Revise: Only use bootstrap param draws where loglike did not go down during numInitSteps
    if(MLEobj$logLik > bestLL) {
      # update the best initial parameter estimates
      tmp = MLEobj$par
      bestinits$Z = tmp$Z
      bestinits$A = tmp$A
      bestinits$R = tmp$R
      bestinits$B = tmp$B
      bestinits$U = tmp$U
      bestinits$Q = tmp$Q
      bestinits$x0 = tmp$x0
      bestinits$V0 = tmp$V0
      bestLL = MLEobj$logLik
    }
   
  } # end numInits loop

  return(bestinits)

}
