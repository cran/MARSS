#######################################################################################################
#   MARSSaic function
#   Adds to marssMLE object:
#   elements in output arg
#   samp.size, num.params
#######################################################################################################
MARSSaic = function(MLEobj, output=c("AIC","AICc"), Options=list(nboot=1000, return.logL.star=FALSE,
  silent=FALSE)) {
  # Options$nboot is the number of bootstrap replicates to do
  # mssm.model is a specified model (needs the structure and parameter elements of the list)
  # output tells what output to produce; this can be a vector  if multiple items should be returned
  #      AIC, AICc, AICbp, AICbb, AICi, boot.params
  # silent is a flag to indicate whether the progress bar should be printed

  if(!is.list(Options)){
    msg = " Options argument must be passed in as a list.\n"
    cat("\nErrors were caught in MARSSaic \n", msg, sep="") 
    stop("Stopped in MARSSaic() due to problem(s) with arguments.\n", call.=FALSE)
  } 
  ## Set options if some not passed in
  if(is.null(Options$nboot)) Options$nboot = 1000
  if(is.null(Options$return.logL.star)) Options$return.logL.star = FALSE
  if(is.null(Options$silent)) Options$silent = FALSE
  tmp = is.marssMLE(MLEobj)
  if(!isTRUE(tmp)) {
    if(!Options$silent) cat(tmp)
    stop("Stopped in MARSSaic() due to problem(s) with the MLE object passed in.\n", call.=FALSE)
    }
    
  if(is.null(MLEobj$logLik)){
    msg = " No log likelihood.  This function expects a model fitted via maximum-likelihood.\n"
    cat("\n","Errors were caught in MARSSaic \n", msg, sep="") 
    stop("Stopped in MARSSaic() due to problem(s) with the MLE object passed in.\n", call.=FALSE)
    }
  return.list=list()

  ## Some renaming for readability
  model = MLEobj$model
  kf = MLEobj$kf
  loglike = MLEobj$logLik
  
  ##### AIC and AICc calculations
  if("AIC" %in% output | "AICc" %in% output){
    K = 0
    for (elem in c("Z", "A", "B", "U", "x0" )) {
      tmp = as.vector(model$free[[elem]])
      pars = length( unique(tmp[!is.na(tmp)]) )
      K = K + pars
    }
    ## Q, R are varcov matrices
    for (elem in c("R", "Q")) {
      tmp = model$free[[elem]]
      tmp1 = as.vector(tmp[upper.tri(tmp, diag=TRUE)])
      pars = length( unique(tmp1[!is.na(tmp1)]) )
      K = K + pars
    }     
    MLEobj$AIC = -2*loglike + 2*K
    samp.size = ifelse(is.na(model$miss.value), sum(!is.na(model$data)), sum(model$data != model$miss.value) )
    MLEobj$AICc = ifelse(samp.size > (K+1),-2*loglike + 2*K*(samp.size/(samp.size-K-1)),"NA, number of data points less than K+1")
    MLEobj$samp.size = samp.size 
    MLEobj$num.params = K
  }

  ##### AICbb & AICbp  
  method = c("AICbb","AICbp") 

  for (m in method[which(method %in% output)]) {
    bootstrap.method = switch(m,
			AICbb = "innovations",
			AICbp = "parametric")
    logL.star=0

    drawProgressBar = FALSE #If the time library is not installed, no prog bar
    if(!Options$silent) { #then we can draw a progress bar
       cat(paste(m, "calculation in progress...\n"))
       prev = progressBar()
       drawProgressBar = TRUE
    }

    boot.params = MARSSboot(MLEobj, nboot=Options$nboot, output="parameters", sim=bootstrap.method,
          param.gen="KalmanEM", silent=TRUE)$boot.params

    # feature request 540
    if((bootstrap.method == "parametric") && ("boot.params" %in% output)) MLEobj$boot.params = boot.params

    for(i in 1:Options$nboot) {
      boot.model = MARSSvectorizeparam( MLEobj, parvec=boot.params[,i] )
      logL.star[i] = MARSSkf(boot.model$model$data, boot.model$par, missing.matrix=boot.model$model$M, miss.value=boot.model$model$miss.value)$logLik	
      if(drawProgressBar) prev = progressBar(i/Options$nboot,prev)
    }
    MLEobj[[m]] = -4*(sum(logL.star))/Options$nboot + 2*loglike   # -2*model$loglik + 2*(1/N)*(-2)*sum(boot$Yloglike-model$logLik)
    if(Options$return.logL.star == TRUE) {
      tmp = paste(m, ".logL.star", sep="")
      MLEobj[[tmp]] = logL.star
    }
  }

  ##### AICbi  using Bengstrom and Cavanaugh algorithm
  if("AICi" %in% output){
    stop("Stopped in MARSSaic() because AICi computation is not fully tested at this time.\n", call.=FALSE)
    #    I wasn't able to replicate Bengstrom and Cavanaugh's results with my code
    logL.AICi=0; logL2.AICi=0; eqn8=0;
    if(!Options$silent) { #then we can draw a progress bar
      cat("AICi calculation in progress...\n")
      prev = progressBar()
      drawProgressBar = TRUE
    }
  TT = dim(model$data)[1]
  n = dim(MLEobj$par$Z)[1]
  m = dim(MLEobj$par$Z)[2]
  for(i in 1:Options$nboot){
      simData = rmvnorm(TT, mean = rep(0, n), sigma = diag(n), method="chol")  #make 0,1
      Astart = rep(0,n)
      Rstart = rep(0.5,n)
      Qstart = rep(0.5,m)
      Ustart = rep(0, m)
      x0start = rep(0, m)
      # estimate the MLEs (theta*) using the 0,1 sim data
      newmod = MLEobj$model
      boot.control = MLEobj$control; boot.control$silent = TRUE  #keep MARSSkem from printing output at each boot
      newmod$data = t(simData)
      start = MLEobj$start
      start$A = Astart; start$R = Rstart; start$Q = Qstart; start$U = Ustart; start$x0 = x0start
      newobj = list(model=newmod, start=start, control=boot.control )
      boot.MLEobj = MARSSkem(newobj) 
      # store the logL of theta* under the sim data
      logL.AICi[i] = boot.MLEobj$logLik
      # generate more sim data
      simData2 = rmvnorm(TT, mean = rep(0, n), sigma = diag(n), method="chol")  #make 0,1
      # obtain the logL of theta* under the sim.star data
      kf.star = MARSSkf(y=t(simData2), parList=boot.MLEobj$par, miss.value=-99) #no missing values
      logL2.AICi[i] = kf.star$logLik

      #B&C eqn 8
      tmpt = 0  # just a temporary holder
      for(t in 1:TT) {
          include.Y <- boot.MLEobj$model$M[,,t]
          dim(include.Y) <- c(n,m)  #stops R from changing the dimensions when m=1
          include.Y2 <- (include.Y%*%array(1,dim=c(n,1))) # make an include.Y for A that is nx1
          tmp = include.Y %*% kf.star$xtt1[,t] + include.Y2*boot.MLEobj$par$A
          Ftinv <- chol2inv(chol(boot.MLEobj$kf$Sigma[,,t]))
          tmpt = tmpt + t(tmp) %*% Ftinv %*% tmp  + sum(diag(Ftinv))
          }
      eqn8[i] = tmpt
      if(drawProgressBar) prev <- progressBar(i/Options$nboot,prev)
      }  #end nboot
      MLEobj$AICi = -2*loglike + (-2)*(1/Options$nboot)*sum(logL2.AICi - logL.AICi)
      Tconstant = TT*n
      MLEobj$AICia = -2 * loglike + (-1 * Tconstant + (1/Options$nboot)*sum( eqn8 ))   #an approximation
      if(Options$return.logL.star == TRUE) {
         MLEobj$AICi.logL.star = logL.AICi
         MLEobj$AICi.logL.star2 = logL2.AICi
         }
    }        #end AICi section

return(MLEobj)
}
