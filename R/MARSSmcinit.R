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
  bounds.tmp = list(B=tmp$B, U=tmp$U, Q=tmp$Q, Z=tmp$Z, A=tmp$A, R=tmp$R)
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
        tmp=as.design(MLEobj$model$fixed[[el]],MLEobj$model$free[[el]])
        if(el %in% c("Q", "R")){   # random starts drawn from a wishart dist
	        if( bounds.param[1] < dim.param[1]){ df=dim.param[1] }else{ df=bounds.param[1] }
	        S=diag(bounds.param[2],dim.param[1])
	        #draw a random matrix from wishart
	        tmp.random = rwishart(df, S)/df
	        #reapply the sharing and fixed constraints 
          element.random = solve(t(tmp$D)%*%tmp$D)%*%t(tmp$D)%*%(vec(tmp.random)-tmp$f)
          param.random = unvec(tmp$f + tmp$D%*%element.random, dim.param)
	      }else{
	        element.random = matrix(runif(dim(tmp$D)[2], bounds.param[1], bounds.param[2]), dim(tmp$D)[2],1)
          param.random = unvec(tmp$f+tmp$D%*%element.random, dim.param)
          if(el %in% c("B")){
           tmp.max=max(abs(eigen(param.random,only.values=TRUE)$values))
           #rescale to bring the max abs eigenvalues to between 0 and 1
           param.random = unvec(tmp$f+tmp$D%*%( element.random/(tmp.max/runif(1,.01,.99)) ), dim.param)
          }
        } 
      }else{ param.random=MLEobj$model$fixed[[el]] }
      init.loop[[el]] = param.random 
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
