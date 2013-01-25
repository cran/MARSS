#######################################################################################################
#   MARSShatyt function
#   Expectations involving hatyt
#######################################################################################################
MARSShatyt = function( MLEobj ) {
    modelObj = MLEobj$model
    if(!is.null(MLEobj[["kf"]])){ kfList = MLEobj$kf
    }else{ kfList=MARSSkf(MLEobj) }
    n=dim(modelObj$data)[1]; TT=dim(modelObj$data)[2]; m=dim(modelObj$fixed$x0)[1]

    #create the YM matrix
    if(is.na(modelObj$miss.value)){ YM=matrix(as.numeric(!is.na(modelObj$data)),n,TT)
      }else{ YM=matrix(as.numeric(!(modelObj$data==modelObj$miss.value)),n,TT) } 
    #Make sure the missing vals in y are zeroed out if there are any
    y=modelObj$data
    y[YM==0]=0

    #set-up matrices for hatxt for 1:TT and 1:TT-1
    IIz=list(); IIz$V0=makediag(as.numeric(takediag(parmat(MLEobj, "V0", t=1)$V0)==0),m)
    hatxt=kfList$xtT      #1:TT
    E.x0 = (diag(1,m)-IIz$V0)%*%kfList$x0T+IIz$V0%*%parmat(MLEobj, "x0", t=1)$x0  #0:T-1
    hatxt1=cbind(E.x0,kfList$xtT[,1:(TT-1),drop=FALSE])
    hatVt=kfList$VtT
    hatVtt1=kfList$Vtt1T

	  msg=NULL

    #Construct needed identity matrices
    I.n = diag(1,n)
	        
    #Note diff in param names from S&S;B=Phi, Z=A, A not in S&S
    time.varying = c(); pari=list()
    for(elem in c("R","Z","A")){  #only params needed for this function
      if( (dim(modelObj$free[[elem]])[3] == 1) & (dim(modelObj$fixed[[elem]])[3] == 1)){  #not time-varying
        pari[[elem]]=parmat(MLEobj, elem, t=1)[[elem]]
        if(elem=="R"){ 
          if(length(pari$R)==1) diag.R=unname(pari$R) else diag.R = takediag(unname(pari$R))
          is.R.diagonal = is.diagonal(pari$R) }
	  }else{ time.varying = c(time.varying, elem) } #which elements are time varying
	  }#end for loop over elem
            
             
    #initialize - these are for the forward, Kalman, filter
    # for notation purposes, 't' represents current point in time, 'TT' represents the length of the series
    hatyt = matrix(0,n,TT)     
    hatOt = array(0,dim=c(n,n,TT))     
    hatyxt = hatyxtt1 = array(0,dim=c(n,m,TT))
    
    for (t in 1:TT) {
      for(elem in time.varying){
        pari[[elem]]=parmat(MLEobj, elem, t=t)[[elem]] 
        if(elem=="R"){ 
          if(length(pari$R)==1) diag.R=unname(pari$R) else diag.R = takediag(unname(pari$R))
          is.R.diagonal = is.diagonal(pari$R) }
      }
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
          Rinv = try(chol(mho.r%*%pari$R%*%t.mho.r))
          #Catch errors before entering chol2inv
          if(class(Rinv)=="try-error") {
            return(list(ok=FALSE, errors="Stopped in MARSShatyt: chol(R) error.\n" ) )      
          }
          Rinv=chol2inv(Rinv) 
          Delta.r = I.n- pari$R%*%t.mho.r%*%Rinv%*%mho.r
        }
        hatyt[,t]=y[,t,drop=FALSE] - Delta.r%*%(y[,t,drop=FALSE]-pari$Z%*%hatxt[,t,drop=FALSE]-pari$A)
        t.DZ = matrix(Delta.r%*%pari$Z,m,n,byrow=TRUE)
        hatOt[,,t]=I.2%*%(Delta.r%*%pari$R+Delta.r%*%pari$Z%*%hatVt[,,t]%*%t.DZ)%*%I.2 + hatyt[,t,drop=FALSE]%*%matrix(hatyt[,t,drop=FALSE],1,n)
        hatyxt[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt[,t,drop=FALSE],1,m)+Delta.r%*%pari$Z%*%hatVt[,,t]
        hatyxtt1[,,t]=hatyt[,t,drop=FALSE]%*%matrix(hatxt1[,t,drop=FALSE],1,m)+Delta.r%*%pari$Z%*%hatVtt1[,,t]
      }
    } #for loop over time
    rtn.list=list(ytT = hatyt, OtT = hatOt, yxtT=hatyxt, yxtt1T=hatyxtt1) 
    return(c(rtn.list,list(ok=TRUE, errors = msg)))
}
