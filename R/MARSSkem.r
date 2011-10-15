#######################################################################################################
#   KEM function
#   Minimal error checking is done.  You should run is.marssMLE(MLEobj) and is.marssm(MLEobj$model) before calling this.
#   Maximization using an EM algorithm with Kalman filter
#######################################################################################################
MARSSkem = function(MLEobj) {
# This is a core function and does not check if user specified a legal or solveable model. 
# free and fixed are a list of model matrices.  Values that are not fixed must be designated NA in 'fixed'; 
# Values that are not free must be designated NA in 'free'
# y is MLEobj$model$data with the missing values replace by 0
kf.x0 = MLEobj$control$kf.x0  #the initial conditions treatment "x00" x0 is at t=0 or "x01" x0 is at t=1
#kf.x0=x00 prior is defined as being E[x(t=0)|y(t=0)]; xtt[0]=x0; Vtt[0]=V0
#kf.x1=x10 prior is defined as being E[x(t=0)|y(t=0)]; xtt1[1]=x0; Vtt1[1]=V0

  #Check that model is allowed given the EM algorithm constaints; returns some info on the model structure
  constr.type=MARSSkemcheck(MLEobj$model, method=MLEobj$method, kf.x0=kf.x0)
  #set up holders for warning messages
  msg=NULL; stop.msg=NULL; msg.kem=NULL; msg.kf=NULL; msg.conv=NULL #error messages

  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MLEobj$model$data #must have time going across columns
  model = MLEobj$model
  free = MLEobj$model$free
  fixed = MLEobj$model$fixed
  inits = MLEobj$start
  model.el = names(fixed)
  n = dim(y)[1]; TT = dim(y)[2]; m = dim(as.matrix(fixed$Q))[1]
  Id = list(m = diag(1,m), n = diag(1,n)) # identity matrices
  control = MLEobj$control
  if(identical(unname(MLEobj$model$fixed$V0), matrix(0,m,m))) x0.is.stochastic = FALSE else x0.is.stochastic = TRUE

  stopped.with.errors=FALSE; kf=NULL; condition.limit=1E10
  tmp=takediag(MLEobj$model$fixed$V0)
  ## Warn about small initial value for V0
  if(any(!is.na(tmp)) && !is.fixed(MLEobj$model$fixed$x0) && x0.is.stochastic && any(tmp[!is.na(tmp)] < 0.05) ){
    msg.conv=c(msg.conv,"One of the diagonal elements of V0 is small. This will cause x0 to converge slowly . See ?MARSSkem.\n")
  }
  tmp = takediag(MLEobj$model$fixed$Q)
  if(any(!is.na(tmp)) && any(tmp[!is.na(tmp)]<0.01) && !is.fixed(MLEobj$model$fixed$U) ){
    msg.conv=c(msg.conv, "Some Q on the diagonal are fixed very small (<0.01). This could cause U to converge slowly.\n")
  }
       
  ## assign the starting parameter values; use fixed values where fixed otherwise use inits; V0 will be reassigned below
  ## do the test for whether parameter is fixed here
  for(elem in model.el) {
    inits[[elem]][is.na(free[[elem]])] = fixed[[elem]][is.na(free[[elem]])]
    assign(elem, inits[[elem]])
    } 
  
  ## create fixed0 matrices with NAs replaced with 0s so we can add the estimated to fixed to get the $par matrices
  fixed0 = fixed
  for(elem in model.el) fixed0[[elem]][is.na(fixed[[elem]])] = 0
  
  ## create the fixec vec and design matrices for each parameter
  d = f = not.fixed = list()
  for(elem in model.el) {
    design = as.design(fixed[[elem]],free[[elem]])
    d[[elem]] = design$D
    f[[elem]] = design$f
    not.fixed[[elem]]=!is.fixed(fixed[[elem]])
    }
      
  ## M is the matrix for handling missing values; 0 on diag = missing
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
  conv.test=list(convergence=72, messages="No convergence test performed.\n", not.converged.params=model.el, converged.params=c() )   # 72 means no info yet; 0 means converged
  cvg = ifelse(is.null(control$abstol), 1, 1 + control$abstol )
  cvg2 = ifelse(is.null(control$abstol), 1, 1 + control$abstol )
  loglike.new = NA #start with no value

  for(iter in 1:control$maxit) { 
    ################# E STEP Estimate states given U,Q,A,R,B,X0 via Kalman filter
    #####################################################################################
    iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )    #the parameter list at iteration iter
    kf.last = kf
    kf = MARSSkf(y, iter.params, missing.matrix=M, init.state=kf.x0, debugkf=control$trace)   #kf.x0 set at top
    if(control$demean.states) {
      xbar = apply(cbind(kf$x0T,kf$xtT),1,mean)
      kf$xtT = kf$xtT-xbar
      kf$x0T = kf$x0T-xbar
    }
    if(!kf$ok) { 
      if(control$trace){ msg.kf=c(msg.kf,paste("iter=",iter," ",kf$errors) )
      }else msg.kf=kf$errors
      stop.msg = paste("Stopped at iter=",iter," in MARSSkem() because numerical errors were generated in MARSSkf\n",sep="")
      stopped.with.errors=TRUE; break
      }
    Ey = MARSShatyt(y, iter.params, kf, missing.matrix=M)
    if(!Ey$ok) { 
      if(control$trace){ msg.kf=c(msg.kf,paste("iter=",iter," ",Ey$errors) )
      }else msg.kf=Ey$errors
      stop.msg = paste("Stopped at iter=",iter," in MARSSkem() because numerical errors were generated in MARSShatyt\n",sep="")
      stopped.with.errors=TRUE; break
      }
    loglike.new = kf$logLik
   
   # This is a diagnostic line that checks if the solution is becoming unstable
    if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg = loglike.new - loglike.old  
    if(iter > 2 & cvg < -sqrt(.Machine$double.eps)) {
      if(control$trace){ 
          msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED.  old=", loglike.old, " new=", loglike.new, "\n", sep=""))
      }else msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED.\n"
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
    # Update Q and R
    # Run Kalman smoother again to update the hidden states expectations
    # Update the other parameters
    
    ################################################################
    # Get new R subject to its constraints
    ################################################################
    #Start the testing for 0s along the diagonal of R
    elem = "R"; thedim = n
    if(not.fixed[[elem]] && control$allow.degen && iter>control$min.degen.iter){ 
      names.diag = paste(elem,free[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)],sep=".") #create a name list that is like that in iter.record
      names.diag[is.na(free[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)])]=NA   #any NAs should be NA not Q.NA
      diag.param = get(elem)[1L + 0L:(thedim - 1L) * (thedim + 1L)]
      diag.fixed = fixed[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)]
      if(any(diag.param<control$degen.lim & is.na(diag.fixed)) && iter>=control$min.iter.conv.test && iter<control$minit) # run conv test 
          conv.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
      degen.elements = diag.param<control$degen.lim & is.na(diag.fixed) & names.diag%in%conv.test$not.converged.param
      if( any( degen.elements ) ){
        #if R=0, then corresponding A and Z rows must be fixed
        allowed.to.be.degen = !is.na(fixed$A) & !apply(is.na(fixed$Z),1,any)
        if(any(!allowed.to.be.degen & degen.elements)) msg.kf=c(msg.kf,paste("Warning: At iter=",iter," attempt to set 0 diagonals for R blocked for elements where corresponding rows of A or Z are not fixed.\n", sep="") ); 
         degen.elements = degen.elements & allowed.to.be.degen
      }
      if( any( degen.elements ) ){
        #if R=0, then corresponding A and Z rows must be fixed
        degen.param=get(elem)
        degen.param[degen.elements,]=0; degen.param[,degen.elements]=0  #set variances and corr covariances to 0         
        current.params =  list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
        current.params[[elem]]=degen.param
        loglike.old = loglike.new
        new.kf = MARSSkf(y, current.params, missing.matrix=M, init.state=kf.x0, debugkf=control$trace)
        if(!new.kf$ok) msg.kf=c(msg.kf,paste("Warning: kf returned error at iter=",iter," in attempt to set 0 diagonals for ", elem,"\n", new.kf$errors,"\n Perhaps Q and R are both going to 0?\n", sep="") ); 
        if(new.kf$ok && is.finite(loglike.old) && is.finite(new.kf$logLik) ) tmp.cvg2 = new.kf$logLik - loglike.old  else tmp.cvg2=Inf
        if(new.kf$ok && tmp.cvg2 < -sqrt(.Machine$double.eps)) {
          msg.kem=c(msg.kem,paste("Warning: setting diagonal to 0 blocked at iter=",iter,". logLik was lower in attempt to set 0 diagonals on ",elem," logLik old=", loglike.old, " new=", new.kf$logLik,"\n", sep=""))
        }
        if(new.kf$ok && tmp.cvg2 > -sqrt(.Machine$double.eps)) { #this means degenerate R has lower LL, so accept it
          assign(elem, degen.param)
          fixed[[elem]][degen.elements,]=0; fixed[[elem]][,degen.elements]=0           
          free[[elem]][degen.elements,]=NA; free[[elem]][,degen.elements]=NA           
          design = as.design(fixed[[elem]],free[[elem]]); d[[elem]] = design$D; f[[elem]] = design$f
          not.fixed[[elem]]=!is.fixed(fixed[[elem]])
          #update the kf and hatyt values with the new Q with 0s         
          kf=new.kf
          if(control$demean.states) {
            xbar = apply(cbind(kf$x0T,kf$xtT),1,mean)
            kf$xtT = kf$xtT-xbar
            kf$x0T = kf$x0T-xbar
          }
          Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
          cvg2=tmp.cvg2
          loglike.new=kf$logLik
        }
      }
    }
    #Now run the standard EM update equations
    R.update = 0
    if(not.fixed$R){
      sum1 = 0
      t.Z = t(Z); t.A = t(A) #pull out of for loop
      t.hatxt = t(kf$xtT); t.hatyt = t(Ey$ytT)
      for (i in 1:TT) {
        hatyt = Ey$ytT[,i,drop=FALSE]; hatyxt=matrix(Ey$yxtT[,,i],n,m); hatOt = Ey$OtT[,,i]
        hatPt = kf$VtT[,,i]+kf$xtT[,i,drop=FALSE]%*%t.hatxt[i,,drop=FALSE]
        hatxt = kf$xtT[,i,drop=FALSE] 
        sum1 = sum1 + hatOt - hatyxt%*%t.Z - Z%*%t(hatyxt)- hatyt%*%t.A - A%*%t.hatyt[i,,drop=FALSE] + 
            Z%*%hatPt%*%t.Z + Z%*%hatxt%*%t.A + A%*%t.hatxt[i,,drop=FALSE]%*%t.Z + A%*%t.A
      }      
      sum1 = (sum1+t(sum1))/2 #enforce symmetry
      R.update = chol2inv(chol(t(d$R)%*%d$R))%*%t(d$R)%*%vec(sum1/TT)
      R.update = d$R%*%R.update   #zeros where no estimate
    }
    R = unvec(f$R + R.update, dim=dim(R))

    #Start~~~~~~~~Error checking   
    if(any(eigen(R,symmetric=TRUE,only.values=TRUE)$values<0)) {
      stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. R update is not positive definite.\n",sep="")
      stopped.with.errors=TRUE;  
      break }      
    if(control$safe && !is.fixed(fixed$R) ){  
      current.params =  list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("R", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
       stopped.with.errors=TRUE
       msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
       break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        loglike.new=kf$logLik
        cvg2=new.kf$cvg2
        msg.kem=c(msg.kem, new.kf$msg.kem)
      }
    }
    ################
    # Get new Q subject to its constraints
    ################################################################
    #Start the testing for 0s along the diagonal of Q
    elem = "Q"; thedim = m
    #Test if elements of Q are going to zero; allow.degen must be TRUE and B must be fixed
    if(not.fixed[[elem]] && !not.fixed$B && control$allow.degen && iter>control$min.degen.iter){ 
      names.diag = paste(elem,free[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)],sep=".") #create a name list that is like that in iter.record
      names.diag[is.na(free[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)])]=NA   #any NAs should be NA not Q.NA
      diag.param = get(elem)[1L + 0L:(thedim - 1L) * (thedim + 1L)]
      diag.fixed = fixed[[elem]][1L + 0L:(thedim - 1L) * (thedim + 1L)]
      if(any(diag.param<control$degen.lim & is.na(diag.fixed)) && iter>=control$min.iter.conv.test && iter<control$minit) # run conv test 
          conv.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
      degen.elements = diag.param<control$degen.lim & is.na(diag.fixed) & names.diag%in%conv.test$not.converged.param
      #need to test that B^(0) is diagonal for any u^(0) or x0^(0) that are estimated
      B.degen.and.is.est.u.x0 = B[(is.na(fixed$U) | is.na(fixed$x0)) & degen.elements, (is.na(fixed$U) | is.na(fixed$x0)) & degen.elements]
      B.degen.off.diag1 = B[(is.na(fixed$U) | is.na(fixed$x0)) & degen.elements, !((is.na(fixed$U) | is.na(fixed$x0)) & degen.elements)]
      B.degen.off.diag2 = B[!((is.na(fixed$U) | is.na(fixed$x0)) & degen.elements), (is.na(fixed$U) | is.na(fixed$x0)) & degen.elements]
      if(length(B.degen.off.diag1)==0)  B.degen.off.diag1=B.degen.off.diag2=0
      if( any( degen.elements ) & is.diagonal(B.degen.and.is.est.u.x0) & all(B.degen.off.diag1==0) & all(B.degen.off.diag2==0) ){  
      #if B^0 is not diagonal then cannot make Q^0 0
        degen.param=get(elem)
        degen.param[degen.elements,]=0; degen.param[,degen.elements]=0           
        current.params =  list( Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
        current.params[[elem]]=degen.param
        loglike.old = loglike.new
        new.kf = MARSSkf(y, current.params, missing.matrix=M, init.state=kf.x0, debugkf=control$trace)
        if(!new.kf$ok) msg.kf=c(msg.kf,paste("Warning: kf returned error at iter=",iter," in attempt to set 0 diagonals for ", elem,"\n", new.kf$errors,"\n Perhaps Q and R are both going to 0?\n", sep="") ); 
        if(new.kf$ok && is.finite(loglike.old) && is.finite(new.kf$logLik) ) tmp.cvg2 = new.kf$logLik - loglike.old  else tmp.cvg2=Inf
        if(new.kf$ok && tmp.cvg2 < -sqrt(.Machine$double.eps)) {
          msg.kem=c(msg.kem,paste("Warning: setting diagonal to 0 blocked at iter=",iter,". logLik was lower in attempt to set 0 diagonals on ",elem," logLik old=", loglike.old, " new=", new.kf$logLik,"\n", sep=""))

        }
        if(new.kf$ok && tmp.cvg2 > -sqrt(.Machine$double.eps)) { #this means degenerate Q has lower LL, so accept it
          assign(elem, degen.param)
          fixed[[elem]][degen.elements,]=0; fixed[[elem]][,degen.elements]=0           
          free[[elem]][degen.elements,]=NA; free[[elem]][,degen.elements]=NA           
          design = as.design(fixed[[elem]],free[[elem]]); d[[elem]] = design$D; f[[elem]] = design$f
          not.fixed[[elem]]=!is.fixed(fixed[[elem]])
          #update the kf and hatyt values with the new Q with 0s         
          kf=new.kf
          if(control$demean.states) {
            xbar = apply(cbind(kf$x0T,kf$xtT),1,mean)
            kf$xtT = kf$xtT-xbar
            kf$x0T = kf$x0T-xbar
          }
          Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
          cvg2=tmp.cvg2
          loglike.new=kf$logLik
        }
      }
    }
    #Then do the regular EM update
    Q.update = 0
    if( not.fixed$Q ){
      # If you treat x0 as at t=1 then 
      # S00 = 0; S11 = 0; S10 = 0; X1 = 0; X0 = 0; TT.numer = TT-1
      # Otherwise if x0 is at t=0 follow Shumway and Stoffer (S&S2006 Eqn 6.67-69)
    	t.U = t(U); t.B = t(B) #pull out of for loop
      if(x0.is.stochastic){  #meaning x0 is the mean and V0 is the var of x[t=t0]
        X0 = kf$x0T 
        S00 = kf$V0T + kf$x0T%*%t(kf$x0T)
        S10 = matrix(kf$Vtt1T[,,1],m,m) + kf$xtT[,1,drop=FALSE]%*%t(kf$x0T);
      }else{       #meaning x0 is x[t=t0] and V0 is 0
        X0 = x0
        S00 = x0%*%t(x0)
        S10 = kf$xtT[,1,drop=FALSE]%*%t(x0) # kf$Vtt1T[,,1] = 0 in this case
      }
      hatxt = kf$xtT[,1,drop=FALSE]
      t.hatxt = t(kf$xtT)
      S11 = matrix(kf$VtT[,,1],m,m) + kf$xtT[,1,drop=FALSE]%*%t.hatxt[1, ,drop=FALSE]
      X1 = kf$xtT[,1,drop=FALSE]; 
      TT.numer = TT
      if(kf.x0=="x10"){
        S00 = 0; S11 = 0; S10 = 0; X1 = 0; X0 = 0; TT.numer = TT-1
        #Gharamani treatment of initial condition; the initial condition specifies x at t=1
      }
      for (i in 2:TT) {
        S00 = S00 + kf$VtT[,,i-1] + kf$xtT[,i-1,drop=FALSE]%*%t.hatxt[i-1, ,drop=FALSE]   #sum 2:T E(xt1T%*%t(xt1T))
        if(m!=1) S00 = (S00 + t(S00))/2  #enforce symmmetry
        S10 = S10 + kf$Vtt1T[,,i] + kf$xtT[,i,drop=FALSE]%*%t.hatxt[i-1, ,drop=FALSE]     #sum 2:T E(xtT%*%t(xt1T))
        S11 = S11 + kf$VtT[,,i] + kf$xtT[,i,drop=FALSE]%*%t.hatxt[i, ,drop=FALSE]         #sum 2:T E(xtT%*%t(xt1T))
        if(m!=1) S11 = (S11 + t(S11))/2  #enforce symmetry
        X0 = X0 + kf$xtT[,i-1,drop=FALSE]                                          #sum 2:T E(xt1T)
        X1 = X1 + kf$xtT[,i,drop=FALSE]                                            #sum 2:T E(xtT)
      }
      Q.update = (S11 - B%*%t(S10) - S10%*%t.B + B%*%S00%*%t.B
          -U%*%t(X1) - X1%*%t.U + U%*%t(B%*%X0) + B%*%X0%*%t.U)/TT.numer + U%*%t.U  #TT.numer is T-1 if x0 is at t=1
      Q.update = vec(Q.update + t(Q.update))/2 #ensure symmetry  
      Q.update = chol2inv(chol(t(d$Q)%*%d$Q))%*%t(d$Q)%*%Q.update #same as vec((T-1)S) in the derivation write-up
      Q.update = d$Q%*%Q.update  #zeros where no estimated value     
    }
    Q = unvec(f$Q + Q.update, dim=dim(Q))

    #Start~~~~~~~~~~~~Error checking
    if(any(eigen(Q,symmetric=TRUE,only.values=TRUE)$values<0)) {
      stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. Q update is not positive definite.\n",sep="")
      stopped.with.errors=TRUE;  
      break }      
    if( control$safe &&  !is.fixed(fixed$Q) ){
      current.params =  list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("Q", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      }
    }  
    #Set up the variance matrices needed for the degenerate case
    OMGp=OMGz=star=list()  #OMGp referes to the Omega^+ matrices; OMGz is Omega^0
    for(elem in c("Q", "R")){
        if(all(takediag(get(elem))==0)){
           OMGp[[elem]] = matrix(0,1,dim(get(elem))[1])
           OMGz[[elem]] = makediag(1,dim(get(elem))[1])
           star[[elem]] = matrix(0,dim(get(elem))[1],dim(get(elem))[1])
        }else{
          OMGp[[elem]]=makediag(1,dim(get(elem))[1])[takediag(get(elem))!=0, , drop=FALSE]
          if(all(takediag(get(elem))!=0)){ OMGz[[elem]]=matrix(0,1,dim(get(elem))[1]) }else{
            OMGz[[elem]]=makediag(1,dim(get(elem))[1])[takediag(get(elem))==0, , drop=FALSE] }
          tmpInv=try(chol2inv(chol(OMGp[[elem]]%*%get(elem)%*%t(OMGp[[elem]]))))
          if(inherits(tmpInv, "try-error")){ 
            stop.msg = paste("Stopped at iter=",iter," in MARSSkem. Inverse of ",elem, ".plus not possible.\n",sep="")
            stopped.with.errors=TRUE;  break}
          star[[elem]] = t(OMGp[[elem]])%*%tmpInv%*%OMGp[[elem]]
          star[[elem]] = (star[[elem]]+t(star[[elem]]))/2     #enforce symmetry
        }
        if( control$trace && !is.fixed(fixed[[elem]]) && kappa(OMGp[[elem]]%*%get(elem)%*%t(OMGp[[elem]]))>condition.limit)
         msg.kem=c(msg.kem,paste("iter=",iter," Unstable estimates because ",elem," is ill-conditioned. C =", kappa(OMGp[[elem]]%*%get(elem)%*%t(OMGp[[elem]])),"\n", sep=""))
    }
    star$Sigma = t(Z)%*%star$R%*%Z+star$Q
    Id$q0 = Id$m-t(OMGp$Q)%*%OMGp$Q
    Id$qplus = t(OMGp$Q)%*%OMGp$Q
    star$R.diam = star$R.sharp = matrix(0,m,m)
    if(any(Id$q0==1)){
      star$B.zero = star$B.zero.1 = Id$q0%*%B%*%Id$q0
      star$B.zero.1[star$B.zero.1==1]=0   #setting the B==1 to 0, deals with (1-b^t)/(1-b) = 1 case when b=1
      star$I.B.Inv = solve(Id$m-star$B.zero.1) #B is not nec symmetric so chol2inv won't work here
      if(kf.x0=="x10"){ 
        B.diam=Id$q0
        B.sharp = matrix(0,m,m) 
      }else{ 
        B.diam = star$B.zero
        B.sharp = (Id$q0 - star$B.zero.1)%*%star$I.B.Inv
        B.sharp[star$B.zero==1]=1 #by definition
        }      
      star$R.diam = B.diam %*% t(Z)%*%star$R%*%Z%*% B.diam
      star$R.sharp = B.sharp %*% t(Z)%*%star$R%*%Z%*% B.sharp
      for(i in 2:TT){
         if(kf.x0=="x10") ii=i-1 else ii=i
         B.diam = star$B.zero^ii
         B.sharp = (Id$q0 - star$B.zero.1^ii)%*%star$I.B.Inv
         B.sharp[star$B.zero==1]=ii #by definition
         star$R.diam = star$R.diam + B.diam %*% t(Z)%*%star$R%*%Z%*% B.diam
         star$R.sharp = star$R.sharp + B.sharp %*% t(Z)%*%star$R%*%Z%*% B.sharp
         }
    }
        
    ################
    # Get new x0 subject to its constraints
    ################################################################
    x0.update = 0
    if(!is.fixed(fixed$x0)){  # some element needs estimating
       if( x0.is.stochastic & !(kf.x0=="x10" & all(R==0)) ){
          x0.update = kf$x0T
          if(constr.type$x0!="unconstrained"){
            V0inv = chol2inv(chol(V0))          
            V0inv = (V0inv+t(V0inv))/2     #enforce symmetry
            x0.update = chol2inv(chol(t(d$x0)%*%V0inv%*%d$x0))%*%t(d$x0)%*%V0inv%*%(kf$x0T-f$x0)
            x0.update = d$x0%*%x0.update
          }
       }
       if( !x0.is.stochastic & !(kf.x0=="x10" & all(R==0)) ){
          if(kf.x0=="x10"){ denom = chol2inv(chol(t(d$x0)%*%(t(Z)%*%star$R%*%Z + star$R.diam + t(B)%*%star$Q%*%B)%*%d$x0))
            }else{ denom = chol2inv(chol(t(d$x0)%*%(star$R.diam + t(B)%*%star$Q%*%B)%*%d$x0)) }
          numer1 = 0
          if(any(Id$q0==1)){
            if(kf.x0=="x10"){ 
              B.diam=Id$q0
              B.sharp = matrix(0,m,m) 
            }else{ 
              B.diam = star$B.zero
              B.sharp = (Id$q0 - star$B.zero.1)%*%star$I.B.Inv
              B.sharp[star$B.zero==1]=1 #by definition
            }
              numer1 = Id$q0%*%B.diam%*%t(Z)%*%star$R%*%(Ey$ytT[,1,drop=FALSE]-Z%*%Id$qplus%*%kf$xtT[,1,drop=FALSE]-Z%*%B.sharp%*%U 
              - Z%*%B.diam%*%f$x0-A) 
              for(i in 2:TT){
                if(kf.x0=="x10") ii=i-1 else ii=i
                B.diam = star$B.zero^ii
                B.sharp =  Id$q0%*%(Id$m-star$B.zero.1^ii)%*%star$I.B.Inv%*%Id$q0
                B.sharp[star$B.zero==1]=ii #by definition
                numer1 = numer1 + Id$q0%*%B.diam%*%t(Z)%*%star$R%*%(Ey$ytT[,i,drop=FALSE]-Z%*%Id$qplus%*%kf$xtT[,i,drop=FALSE]-Z%*%B.sharp%*%U 
              - Z%*%B.diam%*%f$x0-A)
              }
          }
          if(kf.x0=="x10"){ 
            numer2 = Id$qplus%*%t(B)%*%star$Q%*%(kf$xtT[,2,drop=FALSE] - B%*%f$x0 - U) + t(Z)%*%star$R%*%(Ey$ytT[,1,drop=FALSE] -  Z%*%f$x0 - A)
          }else { 
            numer2 = Id$qplus%*%t(B)%*%star$Q%*%(kf$xtT[,1,drop=FALSE] - B%*%f$x0 - U) }
          x0.update = d$x0%*%denom%*%t(d$x0)%*%(numer1 + numer2)
          if( kf.x0=="x10" & any(takediag(R)==0) ){ 
            denom = try( solve(OMGz$R%*%Z%*%t(OMGz$R)) )
            if(inherits(denom, "try-error")){ 
              stop.msg = paste("Stopped at iter=",iter," in MARSSkem. OMGz$R%*%Z%*%t(OMGz$R) is not invertable in x0 update.\n", sep="")
              stopped.with.errors=TRUE;  break }
            x0.update[takediag(R)==0] = denom%*%OMGz$R%*%Ey$ytT[,1,drop=FALSE]
            #element by element multiplicaton to zero out the fixed x0's
            #multiplying d$x0 by a p x 1 matrix of 1s to get a m x 1 column vec with 0s where the fixed x0's are
            x0.update = (d$x0%*%matrix(1,dim(d$x0)[2],1))*x0.update 
          }
       }
       if(kf.x0=="x10" & all(R==0)){
          denom = try(solve(Z))
          if(inherits(denom, "try-error")){ 
            stop.msg = paste("Stopped at iter=",iter," in MARSSkem. Z is not invertable in x0 update.\n", sep="")
            stopped.with.errors=TRUE;  break }
          x0.update = denom%*%Ey$ytT[,1,drop=FALSE]  #this will be m x 1
          #element by element multiplicaton to zero out the fixed x0's
          #multiplying d$x0 by a p x 1 matrix of 1s to get a m x 1 column vec with 0s where the fixed x0's are
          x0.update = (d$x0%*%matrix(1,dim(d$x0)[2],1))*x0.update 
       }
    }
    x0 = f$x0 + x0.update
    #~~~~~~~~Error checking
    if(control$safe & !is.fixed(fixed$x0) ){ 
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("x0", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      }    
    }
    ################
    # Get new V0 subject to its constraints
    ################################################################
    V0.update = 0
    if(!is.fixed(fixed$V0)){  # some element needs estimating (obviously V0!=0)
      V0.update = vec(kf$V0T)
      if(constr.type$V0!="unconstrained"){
        V0.update = chol2inv(chol(t(d$V0)%*%d$V0))%*%t(d$V0)%*%V0.update
        V0.update = d$V0%*%V0.update
        }
    V0 = unvec(f$V0 + V0.update,dim=dim(V0))
    }
    #~~~~~~~~Error checking
    if(any(eigen(V0,symmetric=TRUE,only.values=TRUE)$values<0)) {
      stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. V0 update is not positive definite.\n",sep="")
      stopped.with.errors=TRUE;  
      break } 
    if(control$safe & !is.fixed(fixed$V0) ){ 
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("V0", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      }    
    }
    
    ################
    # Get new A subject to its constraints (update of R will use this)
    ##############################################################
    A.last.iter = A
    A.update=0     
    if( not.fixed$A ){ #if there is anything to update
      sum1 = (Ey$ytT - Z %*% kf$xtT)%*%matrix(1,dim(kf$xtT)[2],1)-TT*f$A
      numer = t(d$A)%*%star$R%*%sum1
      denom = try(chol2inv(chol( t(d$A)%*%star$R%*%d$A ) ))
      if(inherits(denom, "try-error")){
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem. t(d$A)%*%star$R%*%d$A is not invertable.\n", sep="")
        stopped.with.errors=TRUE;  break }
      A.update = d$A%*%(denom%*%numer)/TT  #this will be 0 where A is fixed
    }
    A = f$A + A.update     #f$A is a matrix with 0 for values that will be updated and fixed values otherwise
    #~~~~~~~~Error checking
    if( control$safe & !is.fixed(fixed$A) ){
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("A", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      } 
    }    
    ################
    # Get new U subject to its constraints (update of Q and B will use this)
    ################################################################
    U.update=0
    if( not.fixed$U ){ #if there is anything to update
      if(x0.is.stochastic) E.x0 = kf$x0T else E.x0 = x0
      if(kf.x0=="x10") numer1 = 0 else numer1 = kf$xtT[,1,drop=FALSE] - B%*%E.x0 - f$U
      for (i in 2:TT) numer1 = numer1 + kf$xtT[,i,drop=FALSE] - B%*%kf$xtT[,i-1,drop=FALSE] - f$U
      numer1 = Id$qplus%*%star$Q%*%numer1
      numer2 = 0
      if(any(Id$q0==1)){
        if(kf.x0=="x10") { B.diam=Id$m; B.sharp=matrix(0,m,m) 
        }else{ B.diam = star$B.zero; B.sharp = Id$q0%*%Id$m%*%Id$q0 }
        numer2 = B.sharp%*%t(Z)%*%star$R%*%(Ey$ytT[,1,drop=FALSE] - Z%*%Id$qplus%*%kf$xtT[,1,drop=FALSE] - Z%*%B.diam%*%E.x0 - Z%*%B.sharp%*%f$U-A)
        for (i in 2:TT){
          if(kf.x0=="x10") ii=i-1 else ii = i  #superscript changes depending on x10 or x00
          B.diam = star$B.zero^ii  
          B.sharp =  Id$q0%*%(Id$m-star$B.zero.1^ii)%*%star$I.B.Inv%*%Id$q0
          B.sharp[star$B.zero==1]=ii #by definition
          numer2 = numer2 + B.sharp%*%t(Z)%*%star$R%*%(Ey$ytT[,i,drop=FALSE] - Z%*%Id$qplus%*%kf$xtT[,i,drop=FALSE] - Z%*%B.diam%*%E.x0 - Z%*%B.sharp%*%f$U-A)
          }    
      }     
      numer = numer1 + numer2
      if(kf.x0=="x10") TT.u = TT-1 else TT.u = TT #summation indexing in 2nd sum is 2:TT if kf.x0=x10
      denom = try(chol2inv(chol( t(d$U)%*%(star$R.sharp + TT.u*star$Q)%*%d$U ) ))
      if(inherits(denom, "try-error")){ 
            stop.msg = paste("Stopped at iter=",iter," in MARSSkem. t(d$U)%*%(star$R.sharp + T*star$Q)%*%d$U is not invertable.\n", sep="")
            stopped.with.errors=TRUE;  break }
      U.update = d$U%*%denom%*%t(d$U)%*%numer     #U.update will be 0 when that value is not updated
    }
    U = f$U + U.update     #f$U is a matrix with 0 for values that will be updated and fixed values otherwise
    #~~~~~~~~Error checking  
    if(control$safe & !is.fixed(fixed$U) ){
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("U", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      } 
    }
    ################
    # Get new B subject to its constraints
    ################################################################
    B.update=0 
    if( not.fixed$B ) {
        t.kf.xtT = t(kf$xtT) #more t() out of for loop
        Qinv = chol2inv(chol(Q))
        if(x0.is.stochastic){   #meaning x0 is mean and V0 is var
          hatxtm =  kf$x0T; t.hatxtm = t(kf$x0T); hatVtm = kf$V0T 
        }else{  #meaning x0 is x[t=t0] and V0=0
          hatxtm =  x0; t.hatxtm = t(x0); hatVtm = V0 }
        if(kf.x0=="x00"){  #prior is defined as being E[x(t=0)|y(t=0)]; xtt[0]=x0; Vtt[0]=V0
          hatxt = kf$xtT[,1,drop=FALSE]
          Ptm = hatVtm + hatxtm%*%t.hatxtm
          Pttm = kf$Vtt1T[,,1] + hatxt%*%t.hatxtm
          Sum1 = kronecker(Ptm,Qinv)
          Sum2 = Qinv%*%(Pttm - U%*%t.hatxtm)
        }else{ #prior is defined as being E[x(t=1)|y(t=0)]; xtt1[1]=x0; Vtt1[1]=V0
          Sum1 = Sum2 = 0 #see Ghahramani and Hinton treatment.  Summation starts at t=2
        }        

        for (i in 2:TT) { 
          hatxtm = kf$xtT[,i-1,drop=FALSE]
          t.hatxtm = t.kf.xtT[i-1, ,drop=FALSE]
          hatVtm = kf$VtT[,,i-1] 
          hatxt = kf$xtT[,i,drop=FALSE]
          Ptm = hatVtm + hatxtm%*%t.hatxtm
          Pttm = kf$Vtt1T[,,i] + hatxt%*%t.hatxtm
          Sum1 = Sum1 + kronecker(Ptm,Qinv)
          Sum2 = Sum2 + Qinv%*%(Pttm - U%*%t.hatxtm)
          }
        denom = try( chol2inv(chol(t(d$B)%*%Sum1%*%d$B)) )
        if(inherits(denom, "try-error")){
            stop.msg = paste("Stopped at iter=",iter," in MARSSkem. t(d$B)%*%Qinv%*%d$B is not invertable.\n", sep="")
            stopped.with.errors=TRUE;  break }
        B.update = denom%*%t(d$B)%*%(vec(Sum2)-Sum1%*%f$B)
        B.update = d$B%*%B.update
    }
    B = unvec(f$B + B.update,dim=dim(B))
    #~~~~~~~~Error checking      
    if(control$safe & not.fixed$B ){
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("B", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
    }else{
      kf=new.kf$kf
      Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
      cvg2=new.kf$cvg2
      loglike.new=kf$logLik
      msg.kem=c(msg.kem, new.kf$msg.kem)
      }
    }
    if( control$trace && !is.fixed(fixed$B) ) {
        Ck = kappa(t(d$B)%*%Sum1%*%d$B)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable B estimate because P_{t-1,t-1} is ill-conditioned. C =",round(Ck), "\n", sep=""))
        if(any(abs(eigen(B,only.values=TRUE)$values)>1)) msg.kem=c(msg.kem,paste("iter=",iter," B update is outside the unit circle.", "\n", sep=""))
    }
 
    ################
    # Get new Z subject to its constraints
    ################################################################
    Z.update=0
    if( not.fixed$Z ){
      hatyxt = matrix(Ey$yxtT[,,1], n, m); #funny array call to prevent R from restructuring dims
      Sum1 = Sum2 = 0
      for (i in 1:TT) {
        Pt = kf$VtT[,,i] + kf$xtT[,i,drop=FALSE]%*%t(kf$xtT[,i,drop=FALSE])
        hatyxt = matrix(Ey$yxtT[,,i],n,m); #to prevent R from restructuring dims
        Sum1 = Sum1 + kronecker( Pt, star$R )
        Sum2 = Sum2 + star$R%*%(hatyxt-A%*%t(kf$xtT[,i,drop=FALSE]))
        }
      denom = try(chol2inv(chol(t(d$Z)%*%Sum1%*%d$Z)))
      if(inherits(denom, "try-error")){ 
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem. chol2inv(chol(t(d$Z)%*%Sum1%*%d$Z)) is not invertable.\n", sep="")
        stopped.with.errors=TRUE;  break }
      Z.update =denom%*%t(d$Z)%*%(vec(Sum2)-Sum1%*%f$Z)            
    }
    Z = unvec(f$Z + d$Z%*%Z.update, dim=dim(Z))     
    #Start~~~~~~~~~~~~Error checking
    if(control$safe & !is.fixed(fixed$Z) ){ 
      current.params = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )
      new.kf=rerun.kf("Z", iter, y, current.params, M, control, loglike.new, cvg2, kf.x0)
      if(!new.kf$ok){
        stopped.with.errors=TRUE
        msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
        break
      }else{
        kf=new.kf$kf
        Ey = MARSShatyt(y, current.params, kf, missing.matrix=M)
        cvg2=new.kf$cvg2
        loglike.new=kf$logLik
        msg.kem=c(msg.kem, new.kf$msg.kem)
      }
    } 
    if( control$trace & !is.fixed(fixed$Z) ) {
        Ck = kappa(Sum1)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable Z estimate because P_{t,t} is ill-conditioned. C =",round(Ck), sep=""))
        }
               
    ## Make sure R didn't change the dimensions of any of the matrices- 01/20/09
    if(is.matrix(B)==FALSE) B = as.matrix(B,nrow=m, ncol=m)
    if(is.matrix(Q)==FALSE) Q = as.matrix(Q,nrow=m, ncol=m)
    if(is.matrix(R)==FALSE) R = as.matrix(R,nrow=n, ncol=n)
    if(is.matrix(Z)==FALSE) Z = as.matrix(Z,nrow=n, ncol=m)
    if(is.matrix(U)==FALSE)U=matrix(U,m,1)
    if(is.matrix(A)==FALSE)A=matrix(A,n,1)
    if(is.matrix(x0)==FALSE)x0=matrix(x0,m,1)
    iter.params=list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0 )    #the parameter list at iteration iter    
  }  # end inner iter loop
  
  #prepare the MLEobj to return which has the elements set here
  MLEobj.return = list();    class(MLEobj.return) = "marssMLE"
  MLEobj.return$control=MLEobj$control
  MLEobj.return$start=MLEobj$start
  MLEobj.return$model=MLEobj$model
  MLEobj.return$iter.record = iter.record
  MLEobj.return$numIter = iter
  MLEobj.return$method = "kem"
  
  if(stopped.with.errors){
    if( control$silent==2 ) cat("Stopped due to numerical instability or errors. Print $errors from output for info or set silent=FALSE.\n")      
    #print brief msg.  Full msg printed if silent=F
    msg=c(stop.msg,"par, kf, states, iter, loglike are the last values before the error.\n")
    if(!control$safe) {
        msg=c(msg,"Try control$safe=TRUE which uses a slower but slightly more robust algorithm.\n")
        }
    if(!control$trace) {
        msg=c(msg,"Use control$trace=1 to generate a detailed error report. See user guide for insight.\n")
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
    if(!is.null(control$abstol) || conv.test$convergence==3){  #loglog test has not been run because either abstol used or min.iter.conv not reached
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
  MLEobj.return$logLik = loglike.new

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
  if((!is.null(msg.kem) || !is.null(msg.kf)) && !control$trace){  msg = c(msg,  "\nUse control$trace=1 to generate a more detailed error report.\n") }
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
       #iter.record.par = cbind(iter.record$par,logLik=exp(iter.record$logLik)) #exp because we don't want the log of the log
       iter.record.par = cbind(iter.record$par,logLik=iter.record$logLik) 
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
        #test.loglog=lm(log(test.par)~log(test.len))
        test.loglog=(log(test.par[length(test.par)])-log(test.par[1]))/(log(test.len[length(test.len)])-log(test.len[1]))
        #test.conv[j]=test.loglog$coef[2]
        test.conv[j]=test.loglog
      }
    }   
  }
  if(any(is.na(test.conv))) {
    msg="The log-log degeneracy test produced NAs. Try using control$abstol instead. See user guide.\n"
    return( list(convergence=-2, messages=msg) )
  } 
  if(!is.null(test.conv) && !any(is.na(test.conv)) && any(abs(test.conv)>tol)){
    msg=paste("Warning: the ",names.iter[abs(test.conv)>tol]," parameter value has not converged.\n")
    return( list(convergence=1, messages=msg, not.converged.params=names.iter[abs(test.conv)>tol], converged.params=names.iter[abs(test.conv)<=tol]) ) 
   }else { return( list(convergence=0, messages=NULL, not.converged.params=names.iter[abs(test.conv)>tol], converged.params=names.iter[abs(test.conv)<=tol] ) ) }  #0 means converged successfully
}

rerun.kf = function(elem, iter, y, parList, M, control, loglike.new, cvg2, kf.x0){    #Start~~~~~~~~Error checking
      loglike.old = loglike.new
      msg.kem=NULL 
      kf = MARSSkf(y, parList, missing.matrix=M, init.state=kf.x0,debugkf=control$trace)
      if(control$demean.states) {
        xbar = apply(cbind(kf$x0T,kf$xtT),1,mean)
        kf$xtT = kf$xtT-xbar
        kf$x0T = kf$x0T-xbar
      }
      if(!kf$ok){ 
          msg.kf=c(msg.kf,paste("iter=",iter," ", elem," update ",kf$errors,sep="") ); 
          stop.msg = paste("Stopped at iter=",iter," in MARSSkem after ", elem," update: numerical errors in MARSSkf\n",sep="")
          return(list(ok=FALSE, msg.kf=msg.kf, stop.msg=stop.msg)) }
      loglike.new = kf$logLik
      if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
      if(iter > 2 & cvg2 < -sqrt(.Machine$double.eps)) {
        if(control$trace){ 
          msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED in ",elem," update. logLik old=", loglike.old, " new=", loglike.new,"\n", sep=""))
        }else msg.kem = paste("MARSSkem: The soln became unstable and logLik DROPPED in the",elem, "updates.\n")
        }
    return(list(kf=kf, cvg2=cvg2, loglike.new=loglike.new, msg.kem=msg.kem, ok=TRUE))
}

