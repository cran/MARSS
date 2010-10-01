########################################################################################################################
#   MARSSinnovationsboot.r
#   This bootstraps multivariate time series data using Stoffer and Wall algorithm
#   It creates bootstrap data via sampling from the standardized innovations matrix
#   In the MARSS code, this is referred to as the nonparametric bootstrap.  Strictly speaking, it is not nonparametric.
########################################################################################################################
MARSSinnovationsboot = function(MLEobj, nboot=1000, minIndx=3 ) {
   if(MLEobj$model$miss.value %in% MLEobj$model$data) 
      stop("Stopped in MARSSinnovationsboot() because this algorithm resamples from the innovations and doesn't allow missing values.\n", call.=FALSE)

   # The following need to be present in model: parameter estimates + Sigma, Kt, Innov; the latter 3 are part of $kf
   if(any(is.null(MLEobj$par), is.null(MLEobj$kf), is.null(MLEobj$model$data)))
      stop("Stopped in MARSSinnovationsboot(). The passed in marssMLE object is missing something (par, data, or kf).\n", call.=FALSE)
   
   ## Rename things for code readability
   model=MLEobj$model
   TT = dim(model$data)[2]
   m = dim(MLEobj$par$Z)[2]   
   n = dim(MLEobj$par$Z)[1]
   kf = MLEobj$kf
      
   ##### Set holders for output
   newData = matrix(NA, n, TT)    #store as written for state-space eqn with time across columns
   newStates = matrix(NA, m, TT+1) #store as written for state-space eqn with time across columns  
   boot.data = array(NA, dim=c(dim(model$data), nboot))  
   boot.states = array(NA, dim=c(m, TT, nboot))  

    # Calculate the sqrt of sigma matrices, so they don't have to be computed 5000+ times
    sigma.Sqrt = array(0, dim=c(n, n, TT))
    BKS = array(0, dim=c(m, n, TT)) 
    for(i in 1:TT) {
      std.innovations = stdInnov(kf$Sigma, kf$Innov)  # standardized innovations; time across rows
      eig = eigen(kf$Sigma[,,i])   # calculate inv sqrt(sigma[1])
      sigma.Sqrt[,,i] = eig$vectors %*% makediag(sqrt(eig$values)) %*% t(eig$vectors)
      BKS[,,i] = MLEobj$par$B %*% kf$Kt[,,i] %*% sigma.Sqrt[,,i]   # this is m x n
    }

   for(i in 1:nboot){
      # Then the bootstrap algorithm proceeds
      # Stoffer & Wall suggest not sampling from innovations 1-3 (not much data)
      minIndx = ifelse(TT > 5, minIndx, 1)
      samp = sample(seq(minIndx+1, TT), size = (TT-minIndx), replace = TRUE)
      e = as.matrix(std.innovations)
      e[,1:minIndx] = as.matrix(std.innovations[,1:minIndx])
      e[,(minIndx+1):TT] = as.matrix(std.innovations[,samp])
      # newStates is a[] in the writeup by EH
      newStates[,1] = MLEobj$par$x0
      #newStates[,2] = B %*% x0 + U + B %*% Kt[,,1] %*% sigma.Sqrt[,,1] %*% e[,1]
      for(t in 1:TT) {  
         newStates[,t+1] = MLEobj$par$B %*% newStates[,t] + MLEobj$par$U + BKS[,,t] %*% e[,t]
         newData[,t] = MLEobj$par$Z %*% newStates[,t] + MLEobj$par$A + sigma.Sqrt[,,t] %*% e[,t]
      }

      newStates = newStates[,2:(TT+1)]
      boot.data[,,i] = newData
      boot.states[,,i] = newStates
      # reset newStates to its original dim
      newStates = matrix(NA, m, TT+1)
   }
     return(list(boot.states=boot.states, boot.data=boot.data, model=model, nboot=nboot))
}

######################################################################################################################
#   stdInnov
######################################################################################################################
stdInnov = function(SIGMA, INNOV) {
   # This function added by EW Nov 3, 2008
   # SIGMA is covariance matrix, E are original innovations
   TT = dim(INNOV)[2]
   n = dim(INNOV)[1]
   SI = matrix(0, n, TT)
   for(i in 1:TT) {
      a = SIGMA[,,i]
      a.eig = eigen(a)
      a.sqrt = a.eig$vectors %*% makediag((sqrt(a.eig$values))) %*% solve(a.eig$vectors)
      SI[,i] = chol2inv(chol(a.sqrt)) %*% INNOV[,i]  # eqn 1, p359 S&S
   }
   SI[which(INNOV == 0)]=0   # make sure things that should be 0 are 0
   return(SI)
}
