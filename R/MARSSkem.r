#######################################################################################################
#   KEM function
#   Minimal error checking is done.  You should run is.marssMLE(MLEobj) and is.marssm(MLEobj$model) before calling this.
#   Maximization using an EM algorithm with Kalman filter
#######################################################################################################
MARSSkem = function(MLEobj) {
# This is a core function and does not check if user specified a legal or solveable model. 
# free and fixed are list of constraint matrices.  Values that are not fixed must be designated NA in 'fixed'; 
# Values that are not free must be designated NA in 'free'
# y is MLEobj$model$data with the missing values replace by 0

  #Check that model is allowed given the EM algorithm constaints; returns some info on the model structure
  constr.type=MARSSkemcheck(MLEobj$model, method=MLEobj$method)
  #set up holders for warning messages
  msg=NULL; stop.msg=NULL; msg.kem=NULL; msg.kf=NULL; msg.conv=NULL #error messages
  
  ## Warn about small initial value for V0
  if(MLEobj$control$iter.V0 < 0.5 & !is.fixed(MLEobj$model$fixed$x0) ){
    msg.conv=c(msg.conv,"control$iter.V0=0 is small. This will cause x0 to converge slowly (and not at all if =0). See ?MARSSkem.\n")
  }
  tmp = takediag(MLEobj$model$fixed$Q)
  if(any(!is.na(tmp)) && any(tmp[!is.na(tmp)]<0.01) && !is.fixed(MLEobj$model$fixed$U) ){
    msg.conv=c(msg.conv, "Some Q on the diagonal are fixed very small (<0.01). This could cause U to converge slowly (and not at all if fixed Q is ca 0).\n")
  }
  if(MLEobj$control$iter.V0 == 0){
    msg.kem=c(msg.kem, "control$iter.V0=0. This makes the EM algorithm unstable. See ?MARSSkem.\n")
  }

  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MLEobj$model$data #must have time going across columns
  model = MLEobj$model
  free = MLEobj$model$free
  fixed = MLEobj$model$fixed
  inits = MLEobj$start
  model.el = names(fixed)
  n = dim(y)[1]; TT = dim(y)[2]; m = dim(as.matrix(fixed$Q))[1]
  control = MLEobj$control
  stopped.with.errors=FALSE; kf=NULL; condition.limit=1E10
     
  ## assign the starting parameter values; use fixed values where fixed otherwise use inits; V0 will be reassigned below
  for(elem in model.el) {
    inits[[elem]][is.na(free[[elem]])] = fixed[[elem]][is.na(free[[elem]])]
    assign(elem, inits[[elem]])
    } 

  ## create fixed0 matrices with NAs replaced with 0s so we can add the estimated to fixed to get the $par matrices
  fixed0 = fixed
  for(elem in model.el) fixed0[[elem]][is.na(fixed[[elem]])] = 0
  
  # If V0 is fixed to be zero, then the EM algorithm is run with a diag V0 set large.  At end, the kalman filter is rerun with 
  if(!is.fixed(fixed$x0)){
   if(!identical(unname(fixed$V0), array(0,dim=c(m,m)))){ 
      stop("Stopped in MARSSkem(). If x0 is estimated, V0 must be 0.  See discussion regarding initial conditions in manual.\n",call.=FALSE)
   }else{
        D=as.design(fixed$x0, free$x0)$D  #need this many places
        V0 = control$iter.V0 * D%*%t(D) #if some x0 are shared, they need V0 with 100% correlation
        }
   }else { V0 = fixed$V0 }   #if x0 is fixed, x0 is treated as a prior
  MLEobj$start$V0=V0 #set to whatever V0 is fixed to

  ## M is the matrix for handling missing values
  M = MLEobj$model$M
  #Make sure the missing vals in y are zeroed out
  for(i in 1:dim(y)[2]){ y[!as.logical(takediag(M[,,i])),i]=0 }

  ## Set up variable for debuging and diagnostics
    debugMLEobj=MLEobj
    debugMLEobj$par=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
    iter.record=list(par=NULL,logLik=NULL)  

  ################# The main EM loop which will run until tol reached or max.iter reached
  #######################################################################################   
  
  #set up the convergence flags
  conv.test=list(convergence=72, messages="No convergence test performed.\n")   # 2 means no info yet; 0 means converged
  cvg = ifelse(is.null(control$abstol), 1, 1 + control$abstol )
  cvg2 = ifelse(is.null(control$abstol), 1, 1 + control$abstol )
  loglike.new = NA #start with no value

  for(iter in 1:control$maxit) { 
    ################# E STEP Estimate states given U,Q,A,R,B,X0 via Kalman filter
    #####################################################################################
    iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )    #this list at iteration iter
    kf.last = kf
    kf <- MARSSkf(y, iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
    if(!kf$ok) { 
      if(control$trace) msg.kf=c(msg.kf,paste("iter=",iter," ",kf$errors) )
      else msg.kf=kf$errors
      stop.msg = paste("Stopped at iter=",iter," in MARSSkem() because numerical errors were generated in MARSSkf\n",sep="")
      stopped.with.errors=TRUE; break
      }
    loglike.new = kf$logLik
   
   # This is a diagnostic line that checks if the solution is becoming unstable
    if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg = loglike.new - loglike.old  
    if(iter > 2 & cvg < -sqrt(.Machine$double.eps)) {
        if(control$trace){ msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED.  old=", loglike.old, " new=", loglike.new, "\n", sep=""))
        }else msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }

    ################
    # Keep a record of the iterations for debugging and convergence diagnostics
    ################################################################
    if(control$trace){ # if trace is on, keep the full record over all iterations
      debugMLEobj$par=iter.params
      iter.record$par=rbind(iter.record$par,MARSSvectorizeparam(debugMLEobj))
      iter.record$logLik=c(iter.record$logLik,loglike.new)
      if(!is.null(kf$errors)) {
        msg.kf=c(msg.kf, paste("iter=",iter," ", kf$errors, sep=""))
        }
    }else { #Otherwise keep just last (control$conv.test.deltaT+1) iterations for diagnostics
      debugMLEobj$par=iter.params
      iter.record$par=rbind(iter.record$par,MARSSvectorizeparam(debugMLEobj))
      iter.record$logLik=c(iter.record$logLik,loglike.new)
      tmp.len=dim(iter.record$par)[1]
      if(tmp.len>(control$conv.test.deltaT+1)) {
        iter.record$par=as.matrix(iter.record$par[(tmp.len-control$conv.test.deltaT):tmp.len,,drop=FALSE])
        iter.record$logLik = iter.record$logLik[(tmp.len-control$conv.test.deltaT):tmp.len]
        }
      }
       
    ################
    # Convergence Test
    ################################################################
    if(iter >= control$minit){  # then do convergence testing
     if(!is.null(control$abstol)){
       if( cvg >= 0 && cvg < control$abstol ){
        conv.test$convergence=0  # means converged
        break  # break out of the iterations loop
        }else conv.test$convergence=1  # means NOT converged
      }else {  # else use the log-log convergence test; the default behavior
        if(iter>=control$min.iter.conv.test){ 
          conv.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
          if(conv.test$convergence!=1) break
        }else conv.test$convergence=3  
      } }
      
    # Store loglike for comparison to new one after parameters are updated
    loglike.old = loglike.new
    
    ################# M STEP update U,Q,A,R,B,X0 via ML given x(t) estimate

    ################
    # Get new x0 subject to its constraints
    # S&S Eqn 4.78
    # Note, if x0 is a known prior, x0.update will be all 0s since it is not updated
    ################################################################
    if(FALSE %in% is.na(free$x0)){  # some element needs estimating
      x0 = array(kf$x0T,dim=c(m,1))
      ## impose grouping and constraints (fixed values)
      tmp=table(free$x0, exclude = c(NA, NaN))
      x0.est.levels=names(tmp)
      x0.numGroups <- length(x0.est.levels)
      Zx0 <- matrix(0,m,x0.numGroups)  
      for(i in x0.est.levels) Zx0[which(as.vector(free$x0)==i),which(x0.est.levels==i)] = 1   #as.vector unzips by column
      x0.element.update = array((t(Zx0) %*% x0)/colSums(Zx0), dim=c(x0.numGroups,1))
      x0.update = array(Zx0%*%x0.element.update, dim=c(m,1))
    }
    else x0.update = 0
    x0 = fixed0$x0 + x0.update
    if(control$safe & is.fixed(fixed$x0) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
      if(!kf$ok){ 
         msg.kf=c(msg.kf,paste("iter=",iter," x0 update ",kf$errors,sep="") ); 
         stop.msg = paste("Stopped at iter=",iter," in MARSSkem after x0 update: numerical errors generated in MARSSkf\n",sep="")
         stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in x0 update. logLik old=", loglike.old, " new=", loglike.new, sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }

    ################
    # Get new A subject to its constraints (update of R will use this)
    ##############################################################
    A.last.iter = A     
    if(FALSE %in% is.na(free$A)){ #if there is anything to update
      tmp=table(free$A, exclude=c(NA,NaN))  #free$A is a numeric matrix with NA for those elements that are not updated
      A.est.levels = names(tmp)
      A.numGroups <- length(A.est.levels)
      ZA = matrix(0, n, A.numGroups)   # matrix to allow shared and fixed growth rates; will be 0 where fixed (called F in my write-up)
      for(i in A.est.levels) ZA[which(free$A==i),which(A.est.levels==i)] <- 1 
      Rinv = chol2inv(chol(R))    # this is calculated here because used twice below      
      Rinv = (Rinv+t(Rinv))/2     #enforce symmetry
      sum1 <- 0        
      for (i in 1:TT) {
        A.if.y.missing = (makediag(1,nrow=n)-M[,,i]) %*% A.last.iter  #A.last.iter if missing otherwise 0
        A.if.y.present = M[,,i] %*% (y[,i] - Z %*% kf$xtT[,i]) #put zeros where values are missing
        sum1 <- sum1 + A.if.y.present + A.if.y.missing
      }	 #end for loop over TT
      numer = t(ZA)%*%Rinv%*%sum1
      denom = chol2inv(chol( t(ZA)%*%Rinv%*%ZA ) )
      A.update = ZA%*%(denom%*%numer)/TT   #this will be 0 where A is fixed
    }
    else A.update=0
    A = fixed0$A + A.update     #fixed0$A is a matrix with 0 for values that will be updated and fixed values otherwise
    #Call kf again if safe = TRUE
    if( control$safe & !is.fixed(fixed$A) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10",debugkf=control$trace)
      if(!kf$ok){ 
          msg.kf=c(msg.kf,paste("iter=",iter," A update ",kf$errors,sep="") ); 
          stop.msg = paste("Stopped at iter=",iter," in MARSSkem after A update: numerical errors generated in MARSSkf\n",sep="")
          stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in A update. logLik old=", loglike.old, " new=", loglike.new, sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
      }
    if( control$trace & !is.fixed(fixed$A) ) {
        Ck = kappa(R)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable A estimate because R is ill-conditioned. C =",round(Ck), sep=""))
        }
    ################
    # Get new U subject to its constraints (update of Q and B will use this)
    ################################################################
    # if some state processes share a u, then we need to take the average across processes sharing a u, 
    # taking into account the variance of each process
    if(FALSE %in% is.na(free$U)){ #if there is anything to update
      X1B0 = 0
      for (i in 2:TT) X1B0 = X1B0 + kf$xtT[,i] - B%*%kf$xtT[,i-1]      
      Qinv = chol2inv(chol(Q))    # this is calculated here because used twice below      
      Qinv = (Qinv+t(Qinv))/2     #enforce symmetry
      ## Construct the constraints matrix ZU
      tmp=table(free$U, exclude=c(NA,NaN))  #free$U is a numeric matrix with NA for those elements that are not updated
      U.est.levels = names(tmp)
      U.numGroups <- length(U.est.levels)
      ZU = matrix(0,m,U.numGroups)   # matrix to allow shared growth rates (called F in my write-up)
      for(i in U.est.levels) ZU[which(free$U==i),which(U.est.levels==i)] <- 1 
      numer = t(ZU)%*%Qinv%*%X1B0
      denom = chol2inv(chol( t(ZU)%*%Qinv%*%ZU ) )
      ## The ZU bit is taking the average across values that are shared
      U.update = ZU%*%(denom%*%numer) / (TT-1)     #U.update will be 0 when that value is not updated
    }
    else U.update=0
    U = fixed0$U + U.update     #fixed0$U is a matrix with 0 for values that will be updated and fixed values otherwise
    if(control$safe & !is.fixed(fixed$U) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10",debugkf=control$trace)
      if(!kf$ok){ 
          msg.kf=c(msg.kf,paste("iter=",iter," U update ",kf$errors, sep="") ); 
          stop.msg = paste("Stopped at iter=",iter," in MARSSkem after U update: numerical errors generated in MARSSkf\n",sep="")
          stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in U update. logLik old=", loglike.old, " new=", loglike.new, sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }
    if( control$trace & !is.fixed(fixed$U) ) {
        Ck = kappa(Q)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable U estimate because Q is ill-conditioned. C =",Ck, sep=""))
        }
        
    ################
    # Get new R subject to its constraints
    # S&S 6.72 with addition of grouping and diagonal constraint variants
    ################################################################
    if(FALSE %in% is.na(free$R)){
      sum1 <- 0
      for (i in 1:TT) {
     	## This is the updating equation for R
      ## EH 7.10.08 THIS NEXT LINE IS THE PART THAT IS PREVENTING ESTIMATION OF R COVARIANCES WHEN THERE ARE MISSING VALUES
        R.if.y.missing <- (makediag(1,nrow=n) - M[,,i])%*%R   
        ## This is going to give 0 if have val and R from last iteration if not
        err <- y[,i]- M[,,i]%*%(Z%*%kf$xtT[,i] + A)  #residuals: this will be zero when y[j] is missing]
        R.if.y.present <- err%*%t(err) + (M[,,i]%*%Z)%*%kf$VtT[,,i]%*%t(M[,,i]%*%Z) 
        ## Updated R estimate if have y data
        sum1 <- sum1 + R.if.y.present + R.if.y.missing 
      } #end for loop
      sum1 = (sum1+t(sum1))/2 #enforce symmetry
      R = sum1/TT #this provides the estimate of the R matrix with diagonal and non-diagonal elements
      ## Now add the constraints.  The unfixed values are updated while the fixed values are fixed
      ## This allows shared values on diagonal and shared covariances
      tmp=table(free$R, exclude=c(NA,NaN))  #free$R is a numeric matrix with NA for those elements that are not updated
      R.est.levels = names(tmp)
      R.numGroups <- length(R.est.levels)
      ZR <- matrix(0,n*n,R.numGroups)  # matrix to allow shared measurement errs
      for(i in R.est.levels) ZR[which(as.vector(free$R)==i),which(R.est.levels==i)] = 1 
      R.element.update = array((t(ZR)%*%array(t(R),dim=c(n*n,1)))/colSums(ZR),dim=c(R.numGroups,1))
      R.update = array(ZR%*%R.element.update, dim=c(n,n))  #R.update is 0 if that element is not updated
    }
    else R.update=0
    R = fixed0$R + R.update     #fixed0$R is a matrix with 0 for values that will be updated and fixed values otherwise
    ##### Catch errors
        if(any(eigen(R)$values<0)) {
          stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable and R update is not positive definite.\n",sep="")
          stopped.with.errors=TRUE;  
          break}      
    ##### Use kf call after each update safe
    if(control$safe & !is.fixed(fixed$R) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10",debugkf=control$trace)
      if(!kf$ok){ 
          msg.kf=c(msg.kf,paste("iter=",iter," R update ",kf$errors,sep="") ); 
          stop.msg = paste("Stopped at iter=",iter," in MARSSkem after R update: numerical errors generated in MARSSkf\n",sep="")
          stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in R update. logLik old=", loglike.old, " new=", loglike.new, sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }
    
    	# S10, S00, and S11 calculations required for Q, B and Z updates
      # This is  different than S&S1 Eqn 4.76 (S&S2 Eqn 6.67-69); theirs is based on the likelihood
    	# as written in Eqn 4.69 (6.64) & and Harvey 4.2.19
    	# We don't have a prior on x0 and what you set V0 at affects the Q M-step calculation
    	# to its detriment.  
    	# Instead I'm using the likelihood calculation from Ghahramani and Hinton which is does 
    	# the sum for the Q bit from 2 to T not 1 to T.  This treats the estimated first state as t=1 not t=0.
      # As a result, x0 and V0 drop out of the likelihood. See notes by EH
    	# This seems to work better although we should be able to rewrite this with V0=0 for the case where x0 is treated
    	# as fixed but unknown.
    	# S&S and Harvey would start here
    	# S00 <- kf$V0T + kf$x0T%*%t(kf$x0T)
      # S11 <- kf$VtT[,,1] + (kf$xtT[,1]-U)%*%t(kf$xtT[,1]-U)
      # S10 <- kf$Vtt1T[,,1] + (kf$xtT[,1]-U)%*%t(kf$x0T);
      # I switched on 7/22/08 to this since it seems less sensitive to V0 and finds Q with max L with lower tol setting
      # with S&S it climbs Q at more slowly and cvg hits tol before max is reached
      # Note that because of this difference, the treatment of x0 in the Kalman filter and smoother is different than in S&S
      ################################################################
    S00 = 0; S11 = 0; S10 = 0; X1 = 0; X0 = 0
    for (i in 2:TT) {
      S00 = S00 + (kf$VtT[,,i-1] + kf$xtT[,i-1]%*%t(kf$xtT[,i-1]));   #sum 2:T E(xt1T%*%t(xt1T))
      S10 = S10 + (kf$Vtt1T[,,i] + kf$xtT[,i]%*%t(kf$xtT[,i-1]));     #sum 2:T E(xtT%*%t(xt1T))
      S11 = S11 + (kf$VtT[,,i] + kf$xtT[,i]%*%t(kf$xtT[,i]));         #sum 2:T E(xtT%*%t(xt1T))
      X0 = X0 + kf$xtT[,i-1]                                          #sum 2:T E(xt1T)
      X1 = X1 + kf$xtT[,i]                                            #sum 2:T E(xtT)
    }

    ################
    # Get new B subject to its constraints
    ################################################################
    # 11-2-09 EEH  There are only 3 options for B (currently)
    # B.is.fixed, B.is.unconstrained, B.is.diagonal (means only diagonal is estimated and off-diagonals == 0)
    ## fixed$B without names
    #d = fixed$B
    #rownames(d) = colnames(d) = NULL
                  
    if( FALSE %in% is.na(free$B) ) {
      ok=FALSE
      if(constr.type$B=="unconstrained") {
        B = (S10-U%*%t(X0))%*%chol2inv(chol(S00)); ok=TRUE  #The unconstrained update equation        
        }
      if(constr.type$B=="diagonal and unequal" || constr.type$B=="scalar") {
        B = makediag(takediag(S10-U%*%t(X0))/takediag(S00)); ok=TRUE
        }	 #the unconstrained diagonal B  update eqn; != to diag of above
      if(!ok) stop("MARSSkem: Code bug. B didn't get updated and that shouldn't happen.")
        ## Now add the constraints and grouping.  This is generic code; doesn't depend on B structure
	      B.est.levels = names(table(as.character(free$B)))
        B.numGroups <- length(B.est.levels)
        ZB <- matrix(0, m*m, B.numGroups)  # matrix to allow shared values
        for(i in B.est.levels) ZB[which(as.vector(free$B)==i),which(B.est.levels==i)] <- 1  
        B.element.update = array((t(ZB)%*%array(B,dim=c(m*m,1)))/colSums(ZB),dim=c(B.numGroups,1))
        B.update=array(ZB%*%B.element.update, dim=c(m,m))  #B.update is 0 if that element is not updated
    }
    else B.update=0
    B = fixed0$B + B.update
    
    if(control$safe & !is.fixed(fixed$B) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
      if(!kf$ok){ 
          msg.kf=c(msg.kf,paste("iter=",iter," B update ",kf$errors,sep="") ); 
          stop.msg = paste("Stopped at iter=",iter," in MARSSkem after B update: numerical errors generated in MARSSkf\n",sep="")
          stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in B update. logLik old=", loglike.old, " new=", loglike.new, sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }
    if( control$trace & !is.fixed(fixed$B) ) {
        Ck = kappa(S00)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable B estimate because P_{t-1,t-1} is ill-conditioned. C =",round(Ck), sep=""))
        }
    
    ################
    # Get new Q subject to its constraints
    ################################################################
    if(FALSE %in% is.na(free$Q)){
      Q = (S11 - B%*%t(S10) - S10%*%t(B) + B%*%S00%*%t(B)
          -U%*%t(X1) - X1%*%t(U) + U%*%t(B%*%X0) + B%*%X0%*%t(U))/(TT-1) + U%*%t(U);    
      ## Now add the constraints and grouping.  The unfixed values are updated while the fixed values are fixed
      tmp=table(free$Q, exclude = c(NA, NaN))
      Q.est.levels=names(tmp)
      Q.numGroups <- length(Q.est.levels)
      ZQ <- matrix(0,m*m,Q.numGroups)  # matrix to allow shared q's
      for(i in Q.est.levels) ZQ[which(as.vector(free$Q)==i),which(Q.est.levels==i)] <- 1   #as.vector unzips by column
      Q.element.update = array((t(ZQ)%*%array(t(Q),dim=c(m*m,1)))/colSums(ZQ),dim=c(Q.numGroups,1))
      Q.update = array(ZQ%*%Q.element.update, dim=c(m,m))
    }
    else Q.update = 0
    Q = fixed0$Q + Q.update     #fixed0$Q is a matrix with 0 for values that will be updated and fixed values otherwise
    ##### Catch errors
        if(any(eigen(Q)$values<0)) {
          stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable and Q update is not positive definite.\n",sep="")
          stopped.with.errors=TRUE;  
          break}      
    ##### Use kf call after each update safe
    if(control$safe & !is.fixed(fixed$Q) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
      if(!kf$ok){ 
        msg.kf=c(msg.kf,paste("iter=",iter," Q update ",kf$errors,sep="") ); 
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem after Q update: numerical errors generated in MARSSkf\n",sep="")
        stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in Q update. logLik old=", loglike.old, " new=", loglike.new, "\n", sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }
    
    ################
    # Get new Z subject to its constraints
    ################################################################
    # 2-17-10 EEH  There are only 3 options for Z (currently)
    # Z.is.fixed, Z.is.unconstrained, Z.is.diagonal (means only diagonal is estimated and off-diagonals == 0)
    if(FALSE %in% is.na(free$Z)){
      ok=FALSE
      Z.S11 = S11 + kf$VtT[,,1] + kf$xtT[,1]%*%t(kf$xtT[,1])
      if(constr.type$Z=="unconstrained") {
        sum1 <- 0
        Z.S11.inv = chol2inv(chol(Z.S11))    #sum 1:TT
        for (i in 1:TT) {
          Z.if.y.missing <- (makediag(1,nrow=n) - M[,,i])%*%Z   #zero out rows where y(t) is present
          Z.if.y.present <- ((y[,i]-M[,,i]%*%A)%*%t(kf$xtT[,i]))%*%Z.S11.inv 
          ## Updated Z estimate if have y data
          sum1 <- sum1 + Z.if.y.present + Z.if.y.missing 
        } #end for loop
        Z = array(sum1, dim=c(n,m)) #this provides the estimate of the Z matrix unconstrained
        ok=TRUE
      }
      if(constr.type$Z=="diagonal and unequal") { #not the same as the diagonal of the unconstrained matrix
        sum1 <- 0
        for (i in 1:TT) {
          Z.if.y.missing <- (makediag(1,nrow=n) - M[,,i])%*%Z   #zero out rows where y(t) is present
          Z.if.y.present <- takediag((y[,i]-M[,,i]%*%A)%*%t(kf$xtT[,i]))/Z.S11 
          sum1 <- sum1 + Z.if.y.present + Z.if.y.missing 
        } #end for loop
        Z = array(sum1, dim=c(n,m)) #this provides the estimate of the Z matrix unconstrained
        ok=TRUE        
        }
      if(!ok) stop("MARSSkem: Code bug. Z didn't get updated and that shouldn't happen.")
      ## Now add the constraints.  This is generic code that doesn't depend on the Z constraints
      tmp=table(free$Z, exclude=c(NA,NaN))  #free$Z is a numeric matrix with NA for those elements that are not updated
      Z.est.levels = names(tmp)
      Z.numGroups <- length(Z.est.levels)
      ZZ <- matrix(0,n*m,Z.numGroups)  # matrix to allow shared values
      for(i in Z.est.levels) ZZ[which(as.vector(free$Z)==i),which(Z.est.levels==i)] = 1 
      Z.element.update = array((t(ZZ)%*%array(Z,dim=c(n*m,1)))/colSums(ZZ),dim=c(Z.numGroups,1))
      Z.update = array(ZZ%*%Z.element.update, dim=c(n,m))  #Z.update is 0 if that element is not updated
    }
    else Z.update=0
    Z = fixed0$Z + Z.update     #fixed0$Z is a matrix with 0 for values that will be updated and fixed values otherwise
    if(control$safe & !is.fixed(fixed$Z) ) {
      updated.iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )  
      kf <- MARSSkf(y, updated.iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
      if(!kf$ok){ 
        msg.kf=c(msg.kf,paste("iter=",iter," Z update ",kf$errors,sep="") ); 
        stop.msg = paste("Stopped at iter=",iter,"in MARSSkem after Z update: numerical errors generated in MARSSkf",sep="")
        stopped.with.errors=TRUE;  break}
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < 0) {
        if(control$trace) msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in Z update. logLik old=", loglike.old, " new=", loglike.new, "\n", sep=""))
        msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED at some point in the algorithm.\n"
        }
    }
    if( control$trace & !is.fixed(fixed$Z) ) {
        Ck = kappa(Z.S11)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable Z estimate because P_{t,t} is ill-conditioned. C =",round(Ck), sep=""))
        }
               
    ## Make sure R didn't change the dimensions of any of the matrices- 01/20/09
    if(is.matrix(B)==FALSE) B = as.matrix(B,nrow=m, ncol=m)
    if(is.matrix(Q)==FALSE) Q = as.matrix(Q,nrow=m, ncol=m)
    if(is.matrix(R)==FALSE) R = as.matrix(R,nrow=n, ncol=n)
    if(is.matrix(Z)==FALSE) Z = as.matrix(Z,nrow=n, ncol=m)
    U=array(U,dim=c(m,1))
    A=array(A,dim=c(n,1))
    x0=array(x0,dim=c(m,1))
    
  }  # end inner iter loop
  
  #prepare the MLEobj to return which has the elements set here
  MLEobj.return = list();    class(MLEobj.return) = "marssMLE"
  MLEobj.return$control=MLEobj$control
  MLEobj.return$start=MLEobj$start
  MLEobj.return$model=MLEobj$model
  MLEobj.return$iter.record = iter.record
  MLEobj.return$numIter = iter
  MLEobj.return$method = "kem"
  
  # Checks if the initial condition is treated as fixed but unknown (in which case $fixed$V0 = 0)
  if(!stopped.with.errors && identical(unname(fixed$V0),array(0,dim=c(m,m)))) {
    iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=fixed$V0 )
    kf = MARSSkf(y, iter.params, missing.matrix=M, init.state="x10", debugkf=control$trace)
    if(!kf$ok){ 
    msg.kf=c(msg.kf,paste("Final MARSSkf call with V0=V0: ",kf$errors,sep="") ); 
        stop.msg = paste("Converged successfully but stopped at final MARSSkf call with V0=V0 in MARSSkem.",sep="")
        stopped.with.errors=TRUE} 
    loglike=kf$logLik
  }
  else loglike=loglike.new
  
  if(stopped.with.errors){
    if( control$silent==2 ) cat("Stopped due to numerical instability or errors. Print $errors from output for info or set silent=FALSE.\n")      
    #print brief msg.  Full msg printed if silent=F
    msg=c(stop.msg,"par, kf, states, iter, loglike are the last values before the error.\n")
    if(!control$safe) {
        msg=c(msg,"Try control$safe=TRUE which uses a slower but slightly more robust algorithm.\n")
        }
    if(!control$trace) {
        msg=c(msg,"Use control$trace=TRUE to generate a detailed error report. See manual for insight.\n")
        }
    ## Attach any algorithm errors to the MLEobj
    if(control$trace && !is.null(msg.kem)) msg=c(msg,"\nMARSSkem errors\n",msg.kem)
    if(control$trace && !is.null(msg.kf)) msg=c(msg,"\nMARSSkf errors\n",msg.kf,"\n")    
    MLEobj.return$errors=msg
        
    MLEobj.return$par=iter.params
    MLEobj.return$kf = kf.last
    MLEobj.return$states = kf.last$xtT
    MLEobj.return$convergence = 52
    MLEobj.return$logLik = loglike.new
    return(MLEobj.return)
    }

  ########### Did not stop with errors 
  ## Set the convergence information
  ## Output depends on how it converged and how iterations were determined
    if(!is.null(control$abstol) || conv.test$convergence==3){  #loglog test has not been run
        loglog.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
    }else loglog.test = conv.test
    
    if(conv.test$convergence!=0 && iter==control$maxit){   # stopped because maxit reached
      MLEobj.return$convergence = ifelse(is.null(control$abstol),10,1)
      msg.conv=loglog.test$messages
      if( !control$silent || control$silent==2 )
         cat(paste("Warning! Reached maxit before parameters converged. Maxit was ",control$maxit,".\n",sep=""))
    }else { # either maxit not reached or reached and converged
        MLEobj.return$convergence=72 #this should be reset somewhere below; here for debugging
        if( conv.test$convergence == 0 ){
          MLEobj.return$convergence=ifelse(loglog.test$convergence==0,0,11)   # loglog test passed
          msg.conv=loglog.test$messages
          if( !control$silent || control$silent==2 ) {
            if(iter==control$minit){ 
               if(loglog.test$convergence==0){ cat(paste("Success! algorithm run for ",iter," iterations and parameters converged.\n",sep=""))
               }else cat(paste("algorithm run for ",iter," iterations, abstol reached, but some parameters have not converged.\n",sep=""))
          }else{
            if(loglog.test$convergence==0){ cat(paste("Success! Parameters converged at ",iter," iterations.\n",sep=""))
            }else cat(paste("abstol reached at ",iter," iterations but some parameters have not converged.\n",sep=""))
          }
          if(control$conv.test.slope.tol>0.1) cat(paste("Alert: conv.test.slope.tol is ",control$conv.test.slope.tol,".\nTest with smaller values (<0.1) to ensure convergence.\n",sep="")) 
          }
        }
       if( loglog.test$convergence < 0 ){
           MLEobj.return$convergence = ifelse(is.null(control$abstol),62,12)   #loglog test returned errors
           msg.conv=loglog.test$messages
           if( !control$silent || control$silent==2 ){
              if(!is.null(control$abstol)){
                 cat(paste("abstol reached at ",iter," iterations but log-log convergence test returned errors.\n",sep=""))
              }else cat(paste("Algorithmm stopped at ",iter," iterations because log-log convergence test returned errors.\n",sep=""))
           }
       }
  }
  if(!is.null(msg.conv)) msg=c(msg, "\nConvergence warnings\n", msg.conv)
  ##############################################################
     
  ## Other misc output
  MLEobj.return$par=iter.params
  MLEobj.return$kf = kf
  MLEobj.return$states = kf$xtT
  MLEobj.return$logLik = loglike

  ## Calculate confidence intervals based on state std errors, see caption of Fig 6.3 (p337) Shumway & Stoffer
  if(!is.null(kf$VtT)){
    if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
    if(m > 1) {
      states.se = matrix(0, nrow=m, ncol=TT)
      for(i in 1:TT) states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
    }
    }else  states.se=NULL
  MLEobj.return$states.se = states.se

  if(!is.null(msg.kem)){ msg.kem=c("\nMARSSkem warnings\n", msg.kem); msg=c(msg, msg.kem) }
  if(!is.null(msg.kf)) { msg.kf=c("\nMARSSkf warnings\n", msg.kf); msg=c(msg, msg.kf) }
  if((!is.null(msg.kem) || !is.null(msg.kf)) && !control$trace){  msg = c(msg,  "\nUse control$trace=TRUE to generate a more detailed error report.\n") }
  if((!is.null(msg.kem) || !is.null(msg.kf)) && (!control$silent || control$silent==2) ){
        cat("Alert: Numerical warnings were generated. Print the $errors element of output to see the warnings.\n")
        }
     
  ## Attach any algorithm errors to the MLEobj
  MLEobj.return$errors=msg

  ## Object returned is list(model, start, control, kf, iter.record, numIter, convergence, logLik, states.se, errors)
  return(MLEobj.return)
}

## Run log-log convergence diagnostics
loglog.conv.test = function(iter.record, iter, params.to.test=c("U","x0","R","Q","A","logLik"), deltaT=9, tol=0.5){
  if( !is.list(iter.record) || !all(c("par","logLik") %in% names(iter.record)) || 
    !any(params.to.test %in% c(names(iter.record$par),names(iter.record))) ||
    length(dim(iter.record$par))!=2 || dim(iter.record$par)[1]<=1 || is.null(colnames(iter.record$par)) ){ 
    msg="par list not a proper list (with par and logLik) or too short for conv test or has no column names.\n"
    return( list(convergence=-1, messages=msg) )
  }else {
    if("logLik" %in% params.to.test){
       iter.record.par = cbind(iter.record$par,logLik=exp(iter.record$logLik)) #exp because we don't want the log of the log
       }else iter.record.par=iter.record$par 
    names.iter=colnames(iter.record.par)
    names.sub=strsplit(names.iter,"\\.")
    num.names = length(names.sub)
    p.elems=NULL
    for(j in 1:num.names)p.elems=c(p.elems,names.sub[[j]][1])
    num.varcov = sum( p.elems %in% params.to.test )
    test.conv=rep(0,num.names)
    for( j in 1:num.names ){
     if( p.elems[j] %in% params.to.test ) {
        test.len2=dim(iter.record.par)[1]
      	test.len1=max(1,test.len2-deltaT)
        test.len=(iter-min(test.len2-1, deltaT)):iter 
        test.par = abs(iter.record.par[test.len1:test.len2,j])
        if(any(test.par==0)) test.par = test.par+1   
        test.loglog=lm(log(test.par)~log(test.len))
        test.conv[j]=test.loglog$coef[2]
      }
    }   
  }
  if(any(is.na(test.conv))) {
    msg="The log-log degeneracy test produced NAs. Try using control$abstol instead. See manual.\n"
    return( list(convergence=-2, messages=msg) )
  } 
  if(!is.null(test.conv) && !any(is.na(test.conv)) && any(abs(test.conv)>tol)){
      msg=paste("Warning: the ",names.iter[abs(test.conv)>tol]," parameter value has not converged.\n")
      return( list(convergence=1, messages=msg) ) 
   }else { return( list(convergence=0, messages=NULL ) ) }
}
