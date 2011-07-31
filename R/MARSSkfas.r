#######################################################################################################
#   MARSSkf function based on Koopman and Durbin's KFAS package
#######################################################################################################
MARSSkfas = function(y, parList, missing.matrix = NULL, miss.value= NULL, init.state="x10", debugkf=FALSE, diffuse=FALSE) {
  if(is.null(missing.matrix) && is.null(miss.value)) stop("Stopped in MARSSkfas() because either missing.matrix or miss.value must be specified.\n")
  if(!is.null(miss.value) && !is.na(miss.value) && !is.numeric(miss.value)) stop("Stopped in MARSSkfkd() because miss.value must be numeric (or NA).")
  if(init.state != "x10") stop("Stopped in MARSSkfas() because init.state must be x10. See user guide.\n")
    m=dim(parList$Z)[2]; n=dim(parList$Z)[1]; TT=dim(y)[2]
    a1=rbind(parList$x0,1)
    Zt=cbind(parList$Z,parList$A)
    Ht=parList$R
    Tt=cbind(rbind(parList$B,matrix(0,1,m)),matrix(c(parList$U,1),m+1,1))
    Qt=matrix(0,m+1,m+1); Qt[1:m,1:m]=parList$Q
    Rt=diag(1,m+1)   
    if( diffuse ) { P1inf=matrix(0,m+1,m+1); P1inf[1:m,1:m]=parList$V0; P1=matrix(0,m+1,m+1)
    }else{ P1inf=matrix(0,m+1,m+1); P1=matrix(0,m+1,m+1); P1[1:m,1:m]=parList$V0 }
    
	  M = missing.matrix
    #if the user didn't pass in the M matrix
    #need all the missing values to be speced as NA
    yt=y
    if(is.null(M)){ yt[y==miss.value]=NA 
      }else{ yt[!apply(M,3,diag)]=NA }
    
    kf.out=kf(yt=yt, Zt= Zt, Tt=Tt, Rt=Rt, Ht=Ht, Qt=Qt, a1=a1, P1=P1, P1inf=P1inf, optcal=c(FALSE,FALSE,FALSE,FALSE))
    ks.out=ks(kf.out)
    loglik=-sum(!is.na(yt))/2 * log(2*pi)+kf.out$lik
      
  rtn.list = list(
    xtT = ks.out$ahat,
    VtT = ks.out$Vt, 
    Vtt1T = NULL,
    x0T = ks.out$ahat[,1,drop=FALSE],
    V0T = array(ks.out$Vt[,,1],dim=dim(parList$V0)),
    Vtt = NULL,
    Vtt1 = array(kf.out$Pt[,,1:TT],dim=c(m,m,TT)),
    J=NULL, 
    J0=NULL,
    Kt=NULL, 
    xtt1 = kf.out$at[,1:TT,drop=FALSE], 
    xtt= NULL,
    Innov=NULL, Sigma=NULL,
    logLik=loglik,
    ok=TRUE,
    errors = NULL
    )
}
