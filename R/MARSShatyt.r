#######################################################################################################
#   MARSShatyt function
#   Expectations involving hatyt
#######################################################################################################
MARSShatyt = function(y, parList, kfList, missing.matrix = NULL, miss.value= NULL) {
    if(is.null(missing.matrix) && is.null(miss.value)) stop("Stopped in MARSShatyt() because either missing.matrix or miss.value must be specified.\n")
    if(!is.null(miss.value) && !is.na(miss.value) && !is.numeric(miss.value)) stop("Stopped in MARSShatyt() because miss.value must be numeric (or NA).")

    #Note diff in param names from S&S;B=Phi, Z=A, A not in S&S
    U=parList$U; Q=parList$Q; R=parList$R;
    A=parList$A; x0=parList$x0; V0=parList$V0; B=parList$B;
    Z=parList$Z; #this is the design matrix, called A in S&S  
    n=dim(y)[1]; TT=dim(y)[2]; m=dim(as.matrix(Q))[1]
	  if(length(R)==1) diag.R=unname(R) else diag.R = takediag(unname(R))
	  is.R.diagonal = is.diagonal(R)
    if(identical(unname(V0), matrix(0,m,m))) x0.is.stochastic = FALSE else x0.is.stochastic = TRUE
    hatxt=kfList$xtT
    if(x0.is.stochastic) hatxt1=cbind(kfList$x0T,kfList$xtT[,1:(TT-1),drop=FALSE])
    if(!x0.is.stochastic) hatxt1=cbind(x0,kfList$xtT[,1:(TT-1),drop=FALSE])
    hatVt=kfList$VtT
    hatVtt1=kfList$Vtt1T
	  M = missing.matrix
	  msg=NULL

    #Construct needed identity matrices
    I.n = diag(1,n)
	  
    #if the user didn't pass in the M matrix
    if(is.null(M)){
      if(is.na(miss.value)){ YM=matrix(as.numeric(!is.na(y)),n,TT)
      }else  YM=matrix(as.numeric(!(y==miss.value)),n,TT)
    }else{
      if(n==1){YM=matrix(M,1,TT)
        }else{
          YM = matrix(1,n,TT)
          for(t in 1:TT) YM[,t]=diag(M[,,t]) }
      }
       
    #Make sure the missing vals in y are zeroed out if there are any
    y[YM==0]=0
             
    #initialize - these are for the forward, Kalman, filter
    # for notation purposes, 't' represents current point in time, 'T' represents the length of the series
    hatyt = matrix(0,n,TT)     
    hatOt = array(0,dim=c(n,n,TT))     
    hatyxt = hatyxtt1 = array(0,dim=c(n,m,TT))
    
    for (t in 1:TT) {
      if(all(YM[,t]==1)){  #none missing
          hatyt[,t]=y[,t,drop=FALSE]
          hatOt[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatyt[,t,drop=FALSE],1,n) #matrix() is faster than t()
          hatyxt[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt[,t,drop=FALSE],1,m)
          hatyxtt1[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt1[,t,drop=FALSE],1,m)
      }else{
        I.2 = I.r = I.n; 
        I.2[YM[,t]==1,]=0 #1 if YM=0 and 0 if YM=1
        I.r[YM[,t]==0 | diag.R==0,]=0   #if Y missing or R = 0, then 0
        Delta.r=I.n
        if(is.R.diagonal) Delta.r = I.n-I.r
        if(!is.R.diagonal && any(YM[,t]==1 & diag.R!=0)){           
          mho.r = I.r[YM[,t]==1 & diag.R!=0,,drop=FALSE]
          t.mho.r = I.r[,YM[,t]==1 & diag.R!=0,drop=FALSE]
          Rinv = try(chol(mho.r%*%R%*%t.mho.r))
          #Catch errors before entering chol2inv
          if(class(Rinv)=="try-error") {
            return(list(ok=FALSE, errors="Stopped in MARSShatyt: chol(R) error.\n" ) )      
          }
          Rinv=chol2inv(Rinv) 
          Delta.r = I.n- R%*%t.mho.r%*%Rinv%*%mho.r
        }
        hatyt[,t]=y[,t,drop=FALSE] - Delta.r%*%(y[,t,drop=FALSE]-Z%*%hatxt[,t,drop=FALSE]-A)
        t.DZ = matrix(Delta.r%*%Z,m,n,byrow=TRUE)
        hatOt[,,t]=I.2%*%(Delta.r%*%R+Delta.r%*%Z%*%hatVt[,,t]%*%t.DZ)%*%I.2 + hatyt[,t,drop=FALSE]%*%matrix(hatyt[,t,drop=FALSE],1,n)
        hatyxt[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt[,t,drop=FALSE],1,m)+Delta.r%*%Z%*%hatVt[,,t]
        hatyxtt1[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt1[,t,drop=FALSE],1,m)+Delta.r%*%Z%*%hatVtt1[,,t]
      }
    }
    rtn.list=list(ytT = hatyt, OtT = hatOt, yxtT=hatyxt, yxtt1T=hatyxtt1) 
    return(c(rtn.list,list(ok=TRUE, errors = msg)))
}
