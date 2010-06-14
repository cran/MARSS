## MARSSinits  
## Set up inits
## These will be checked by the MLE object checker.

MARSSinits <- function(modelObj, inits=list(B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=10))
{
if(!is.list(inits)) stop("MARSSinits: inits must be a list")
default = list(B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=10)
for(elem in names(default)){
  if(is.null(inits[[elem]])) inits[[elem]]=default[[elem]]
}
  y = modelObj$data
  m = dim(modelObj$fixed$Z)[2]
  n = dim(modelObj$fixed$Z)[1]
  miss.value = modelObj$miss.value

  D=as.design(modelObj$fixed$x0, modelObj$free$x0)$D  #need this many places
  if(identical(unname(inits$x0),-99) & !is.fixed(modelObj$fixed$x0) ) {
      if(!is.design(modelObj$fixed$Z))
        stop("MARSSinits: because Z is not a design matrix, you'll need to specify inits$x0.  See help file.")
    tmp.Z.mat=modelObj$fixed$Z
    inits$x0 = array(0,dim=c(m,1))
    # do a simple linear regression, estimating the predicted value using the y that x will be scaled to (in case multiple y per x)
    for(i in 1:m) {
      this.index = min(which(tmp.Z.mat[,i]==1))
      y.obs = y[this.index,]
      y.obs[which(y.obs==miss.value)] = NA
      inits$x0[i,1] = lm(y.obs~seq(1:dim(y)[2]))$fitted.values[1]
    }
    #use average if there are shared values
    inits$x0 = D%*%solve(t(D)%*%D)%*%t(D)%*%inits$x0
  }
  if( is.fixed(modelObj$fixed$x0) ){inits$x0 = modelObj$fixed$x0
  }else inits$x0=array(inits$x0,dim=dim(modelObj$fixed$x0))
  
  if(length(inits$V0)==1){
    if( is.fixed(modelObj$fixed$x0) ){ inits$V0=makediag(inits$V0, nrow=m)
      }else inits$V0 = inits$V0 * D%*%t(D) #if some x0 are shared, they need V0 with 100% correlation
    }
  if(length(inits$V0)==m){
    if( is.fixed(modelObj$fixed$x0) || all(inits$V0==0) ) inits$V0=makediag(inits$V0, nrow=m)
    }      
      
  if(!is.matrix(inits$V0)  || dim(inits$V0)[1]!=m || dim(inits$V0)[2]!=m) #there was a problem
     stop("MARSSinits: inits$V0 is not a mxm matrix.  See help file.")
  if(!is.fixed(modelObj$fixed$x0) && !all( (!(D%*%t(D)))*inits$V0==0 ))
     warning("The initial V0 looks wrong (look at start$V0). A wrong init V0 will mean a wrong logLik value.")
  for(elem in c("A","U","Z")) {
    if(is.fixed(modelObj$fixed[[elem]])){ inits[[elem]] = modelObj$fixed[[elem]]
    }else { #use inits but replace any fixed elements with their fixed values
      inits[[elem]]=array(inits[[elem]],dim=dim(modelObj$fixed[[elem]]))
      inits[[elem]][!is.na(modelObj$fixed[[elem]])]= modelObj$fixed[[elem]][!is.na(modelObj$fixed[[elem]])]
    }
    }
  for(elem in c("Q","R","B")) {
    if(is.fixed(modelObj$fixed[[elem]])){ inits[[elem]] = modelObj$fixed[[elem]]
    }else{ if(length(inits[[elem]])==length(modelObj$fixed[[elem]]) ){
      inits[[elem]]= inits[[elem]]=array(inits[[elem]],dim=dim(modelObj$fixed[[elem]]))
    }else{ if(length(inits[[elem]])==1 || length(inits[[elem]])==dim(modelObj$fixed[[elem]])[1]){
      inits[[elem]]=makediag( inits[[elem]], nrow=dim(modelObj$fixed[[elem]])[1] )
      } }
    # replace any fixed elements with their fixed values
    inits[[elem]][!is.na(modelObj$fixed[[elem]])]= modelObj$fixed[[elem]][!is.na(modelObj$fixed[[elem]])]
    }
    if(!is.matrix(inits[[elem]])) #there was a problem.  
     stop(paste("MARSSinits: inits$",elem," is not a matrix.  See help file.",sep=""))
    }
         
  if(is.null(inits$Z) || inits$Z == -99) {
    if(!is.fixed(modelObj$fixed$Z)) {
      stop("MARSSinits: because you are estimating parts of Z, you'll need to specify inits$Z.  See help file.")
    }
    else inits$Z = modelObj$fixed$Z #Z is fixed
  }
  #else user must have passed in inits$Z and it will be checked by MLE object checker      

  inits
}
