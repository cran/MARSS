#######################################################################################################
#   MARSSkf function
#   Kalman filter and smoother
#   ** All eqn refs are to 2nd ed of Shumway & Stoffer (2006): Time Series Analysis and Its Applications
#######################################################################################################
MARSSkf = function(y, parList, missing.matrix ="not given", miss.value="not given", init.state="x10", debugkf=FALSE) {
    if(identical(missing.matrix,"not given") & identical(miss.value,"not given")) stop("MARSSkf: either missing.matrix or miss.value must be specified.")
    if(!(init.state %in% c("x10","x00"))) stop("MARSSkf: init.state must be either x10 or x00. See manual.")
    condition.limit=1E10
    #Set up param names; Phi is used instead of B since S&S use Phi
    U=parList$U; Q=parList$Q; R=parList$R; 
    A=parList$A; Phi=parList$B; x0=parList$x0; V0=parList$V0
    Z=parList$Z; #this is the design matrix, called A in S&S  
    n=dim(y)[1]; TT=dim(y)[2]; m=dim(as.matrix(Q))[1]
	  M = missing.matrix
	  msg=NULL
	  
    #if the user didn't pass in the missing values matrix, construct M from miss.value 
    if(identical(M, "not given")) {
      if(TRUE %in% (y == miss.value)) {
        M <- array(0, dim=c(n,n,TT))  
        for(i in 1:TT)  M[,,i] <- makediag(ifelse(y[,i]!=miss.value,1,0),nrow=n)  
      }
      else  M = array(makediag(1,nrow=n),dim=c(n, n, TT))  
    }
    #Make sure the missing vals in y are zeroed out
    for(i in 1:dim(y)[2]) y[,i]=M[,,i]%*%y[,i]  
         
    #initialize - these are for the forward, Kalman, filter
    # for notation purposes, 't' represents current point in time, 'T' represents the length of the series
    Vtt <- array(0,dim=c(m,m,TT))     # Analagous to S&S Ptt, var[xt,xt|y(1:t)]
    Vtt1 <- array(0,dim=c(m,m,TT))    # Analagous to S&S Ptt1, cov[xt,xt1|y(1:t)]
    xtt <- array(0,dim=c(m,TT))       # E[x(t) | Y(t)]
    xtt1 <- array(0,dim=c(m,TT))      # E[x(t) | Y(t-1)]
    innov <- array(0,dim=c(n,TT))     # these are innovations, parentheses in 6.21 
    vt <- array(0,dim=c(n,TT))        # used for likelihood, vt equivalent to epsilon, eqn 6.62
    Ft <- array(0,dim=c(n,n,TT))      # used for likelihood, Ft equivalent to sigma matrix eqn 6.62
    # these are for the backwards, smoother
    VtT <- array(0,dim=c(m,m,TT))     # var[xt,xt|y(1:T)]
    J <- array(0,dim=c(m,m,TT))       # see eqn 6.49
    Vtt1T <- array(0,dim=c(m,m,TT))   # cov[xt,xt1|y(1:T)]
    xtT <- array(0,dim=c(m,TT))       # E[x | y(1:T)]
    Kt <- array(0, dim=c(m,n,TT))     # 3D matrix of Kalman gain, EW added 11/14/08
    
    ##############################################
    #FORWARD PASS (K filter) gets you E[x(t) given y(1:t)]
    ##############################################
    # In the following, innov and Ft are needed for the likelihood calculation
		# the missing values will contribute 0.0 for the LL calc
    # R_mod is needed for the corrected likelihood calculation when there are missing values
		# See section 12.3 in Time Series: Theory and Methods (1991) Peter J. Brockwell, Richard A. Davis
		# put 1's on the diagonal where there are missing values and zero out the rows and columns

    for (t in 1:TT) {
    #t=1 treatment depends on how you define the initial condition.  Either as x at t=1 or x at t=0
      if(t==1) {
        if(init.state=="x00") {
          xtt1[,1] <- Phi%*%x0 + U   #Shumway and Stoffer treatment of initial states # eqn 6.19   (pi is defined as t=0)
           Vtt1[,,1] <- Phi%*%V0%*%t(Phi) + Q          # eqn 6.20
        }
        if(init.state=="x10") {    #Ghahramani treatment of initial states uses x10 and has no x00 (pi is defined as t=1)
         xtt1[,1] <- x0         
         Vtt1[,,1] <- V0
        }
      }
      else {   #t!=1
       xtt1[,t] <- Phi%*%xtt[,t-1] + U  #xtt1 denotes x_t^(t-1), eqn 6.19
       Vtt1[,,t] <- Phi%*%Vtt[,,t-1]%*%t(Phi) + Q                  # eqn 6.20
      }              
      Vtt1[,,t] <- (Vtt1[,,t]+t(Vtt1[,,t]))/2   #in general Vtt1 is not symmetric but here it is since Vtt and Q are
      siginv1 = M[,,t]%*%Z%*%Vtt1[,,t]%*%t(M[,,t]%*%Z)+R    # bracketed piece of eqn 6.23        
      ####### Catch errors before entering chol2inv
        if(any(takediag(Vtt1[,,t])==0) ) {  #0s on diag of Vtt1 will break the K smoother (below) if t>1
          if(init.state=="x00" || (init.state=="x10" && t>1)){ return(list(ok=FALSE, errors=paste("Stopped in MARSSkf: soln became unstable when zeros appeared on the diagonal of Vtt or Vtt1 at t>1.\n") ) )
          }else if(any(takediag(siginv1)==0)) return(list(ok=FALSE, errors=paste("Stopped in MARSSkf: soln became unstable when zeros appeared on the diagonal of siginv[,,1].\n") ) )
          }
      siginv2=try(chol(siginv1), silent = TRUE)     # now siginv is sig[[i]]^{-1} 
      if(class(siginv2)=="try-error") {
          Ck1 = try(kappa(siginv1))
          Ck1 = ifelse(class(Ck1)=="try-error","Inf",round(Ck1))
          Ck4 = try(kappa(R)) 
          Ck4 = ifelse(class(Ck4)=="try-error","Inf",round(Ck4))
          msg1=paste("Condition num. of siginv[t=",t,"] = ",Ck1," ",sep="")
          msg2=paste("Condition num. of R = ",Ck4," ",sep="")
          return(list(ok=FALSE, errors=paste("Stopped in MARSSkf: chol(Z%*%Vtt1[,,t]%*%t(Z)+R) error. ",msg1,msg2,"\n") ) )       
          }
      ####### Error-checking
      siginv=chol2inv(siginv2)     # now siginv is sig[[i]]^{-1} 
      siginv = (t(siginv)+siginv)/2     #Vtt1 happens to be symmetric since it is V0+Q; although in general E(xt t(xt1)) is not symmetric
      Kt[,,t] <- Vtt1[,,t]%*%t(M[,,t]%*%Z) %*% siginv;    #broke siginv to impose symmetry, eqn 6.23
      innov[,t] <- y[,t] - M[,,t]%*%(Z%*%xtt1[,t] + A)
      xtt[,t] <- xtt1[,t] + Kt[,,t]%*%innov[,t]   # eqn 6.21
      Vtt[,,t] <- Vtt1[,,t]-Kt[,,t]%*%M[,,t]%*%Z%*%Vtt1[,,t]  # eqn 6.22, detail after 6.28
      Vtt[,,t] <- (Vtt[,,t]+t(Vtt[,,t]))/2 #to ensure its symetric
      # Variables needed for the likelihood calculation; see comments above
      R_mod = (diag(n)-M[,,t]) + M[,,t]%*%R%*%M[,,t]
      vt[,t] <- y[,t]- M[,,t]%*%(Z%*%xtt1[,t]+A) #need to hold on to this for loglike calc
      Ft[,,t] <- (M[,,t]%*%Z)%*%Vtt1[,,t]%*%t(M[,,t]%*%Z)+R_mod #need to hold on to this for loglike calc
      Ft[,,t] <- (Ft[,,t]+t(Ft[,,t]))/2 #to ensure its symetric
       
      ####### Error-checking
      if(debugkf) {
          Ck1 = kappa(siginv)
          Ck2 = kappa(Vtt1[,,t])
          Ck3 = kappa(Ft[,,t])
          if(Ck1>condition.limit && !all(Kt[,,t]==0) ) 
          msg=rbind(msg,paste("MARSSkf: solution is becoming unstable.  Condition num. of siginv[t=",t,"] = ",round(Ck1),"\n",sep=""))
          if(Ck2>condition.limit && t>1) 
          msg=rbind(msg,paste("MARSSkf: solution is becoming unstable.  Condition num. of Vtt1[t=",t,"] = ",round(Ck2),"\n",sep=""))
          if(Ck3>condition.limit){
             Ck4 = kappa(R) 
             msg=rbind(msg,paste("MARSSkf: logLik computation is becoming unstable.  Condition num. of Sigma[t=",t,"] = ",round(Ck3)," and of R = ",round(Ck4),".\n",sep=""))
             }
          }
          #Abandon if solution is so unstable that Vtt diagonal became negative
      if(any(takediag(Vtt[,,t])<0) || any(takediag(Vtt1[,,t])<0) ) 
          return(list(ok=FALSE, 
          errors=paste("Stopped in MARSSkf: soln became unstable and negative values appeared on the diagonal of Vtt or Vtt1.\n") ) )
    } #for i to 1:TT
    KT <- Kt[,,t];

    ######################################################
    #BACKWARD PASS (Kalman smoother) gets you E[x(t)|y(1:T)] from E[x(t)|y(1:t)]
    #indexing is 0 to T for the backwards smoother recursions
    xtT[,TT] <- xtt[,TT]  
    VtT[,,TT] <- Vtt[,,TT]
    s <- seq(TT,2)
    for(i in 1:(TT-1)) {
        yr <- s[i]                 #equivalent to T:-1:0  
        Vinv <- chol2inv(chol(Vtt1[,,yr]))
        Vinv <- (Vinv + t(Vinv))/2 #to enforce symmetry after chol2inv call
        J[,,yr-1] <- Vtt[,,yr-1]%*%t(Phi)%*%Vinv     # eqn 6.49
        xtT[,yr-1] <- xtt[,yr-1] + J[,,yr-1]%*%(xtT[,yr]-xtt1[,yr])     # eqn 6.47
        VtT[,,yr-1] <- Vtt[,,yr-1] + J[,,yr-1]%*%(VtT[,,yr]-Vtt1[,,yr])%*%t(J[,,yr-1])  # eqn 6.48
        VtT[,,yr-1] <- (VtT[,,yr-1]+t(VtT[,,yr-1]))/2     #VtT is symmetric
    }

    if(init.state=="x00") { #Shumway and Stoffer treatment of initial conditions
      Vinv <- chol2inv(chol(Vtt1[,,1]))
      Vinv <- (Vinv + t(Vinv))/2 #to enforce symmetry after chol2inv call
      J0 <- V0%*%t(Phi)%*%Vinv                      # eqn 6.49
      x0T <- x0 + J0%*%(xtT[,1]-xtt1[,1]);          # eqn 6.47
      V0T <- V0 + J0%*%(VtT[,,1]-Vtt1[,,1])*t(J0)   # eqn 6.48
      V0T <- (V0T+t(V0T))/2;
    }
    if(init.state=="x10") { #Ghahramani treatment of initial states
      J0 <- J[,,1]
      x0T <- xtT[,1]
      V0T <- VtT[,,1]
    }
    #run another backward recursion to get E[x(t)x(t-1)|y(T)]
    Vtt1T[,,TT] <- (makediag(1,m) - KT%*%M[,,TT]%*%Z)%*%Phi%*%Vtt[,,TT-1] #eqn. 6.55 this is Var(x(T)x(T-1)|y(T))
    s <- seq(TT,3)
    for (i in 1:(TT-2)) {
       yr <- s[i]
       Vtt1T[,,yr-1] <- Vtt[,,yr-1]%*%t(J[,,yr-2]) + J[,,yr-1]%*%(Vtt1T[,,yr]-Phi%*%Vtt[,,yr-1])%*%t(J[,,yr-2])   #eqn 6.56
    }
    if(init.state=="x00") Vtt1T[,,1] <- Vtt[,,1]%*%t(J0) + J[,,1]%*%(Vtt1T[,,2]-Phi%*%Vtt[,,1])%*%t(J0)
    if(init.state=="x10") Vtt1T[,,1] <- NA

    ###########################################################
    #Calculate log likelihood, see eqn 6.62
    #Innovations form of the likelihood
    loglike <- -sum(M)/2*log(2*pi)    #sum(M) is the number of data points
    for (i in 1:TT) {
        if(length(Ft[,,i])==1) detFt <- Ft[,,i] else detFt <- det(Ft[,,i])
        if( detFt<0 || !is.finite(log(detFt)) ){
          return(list(ok=FALSE, Sigma=Ft, errors=paste("Stopped in MARSSkf: log(det(Ft[,,",i,"]))=NA.\n",sep="") ) )
          }
        Ftinv <- chol2inv(chol(Ft[,,i]))
        Ftinv <- (Ftinv +t(Ftinv))/2 #enforce symmetry; Ft is symmetric
        loglike <- loglike - (1/2)%*%t(vt[,i]) %*% Ftinv %*% vt[,i] - (1/2)*log(detFt);
    }
    if( !is.finite(loglike) ) return(list(ok=FALSE, errors=paste("Stopped in MARSSkf: loglike computed to NA.\n") ) )

    return(list(xtT = xtT, VtT = VtT, Vtt1T = Vtt1T, x0T = x0T, V0T = V0T, logLik = loglike, Vtt = Vtt, Vtt1 = Vtt1, J=J, J0=J0, Kt=Kt, xtt1 = xtt1, xtt=xtt, Innov=innov, Sigma=Ft, ok=TRUE, errors = msg))
}
