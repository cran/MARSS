## MARSSinits  
## Set up inits
## These will be checked by the MLE object checker.

MARSSinits <- function(modelObj, inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=5), method)
{
if(!is.list(inits)) stop("Stopped in MARSSinits() because inits must be a list.\n", call.=FALSE)
default = alldefaults[[method]]
for(elem in names(default)){
  if(is.null(inits[[elem]])) inits[[elem]]=default[[elem]]
}
  y = modelObj$data
  m = dim(modelObj$fixed$Z)[2]
  n = dim(modelObj$fixed$Z)[1]
  miss.value = modelObj$miss.value

  D=as.design(modelObj$fixed$x0, modelObj$free$x0)$D  #need this many places
  f=as.design(modelObj$fixed$x0, modelObj$free$x0)$f  #need this many places
  
  if(length(inits$V0)==1){
    if( is.fixed(modelObj$fixed$x0) ){ inits$V0=makediag(inits$V0, nrow=m)
      }else inits$V0 = inits$V0 * D%*%t(D) #if some x0 are shared, they need V0 with 100% correlation
    }
  if(length(inits$V0)==m){
    if( is.fixed(modelObj$fixed$x0) || all(inits$V0==0) ) inits$V0=makediag(inits$V0, nrow=m)
    }      
      
  if(!is.matrix(inits$V0)  || dim(inits$V0)[1]!=m || dim(inits$V0)[2]!=m) #there was a problem
     stop("Stopped in MARSSinits() because inits$V0 is not a mxm matrix.  See help file.\n", call.=FALSE)
  if(!is.fixed(modelObj$fixed$x0) && !all( (!(D%*%t(D)))*inits$V0==0 ))
     warning("The initial V0 looks wrong (look at start$V0). A wrong init V0 will mean a wrong logLik value.")

  for(elem in c("A","U","Z")) {
    if(is.fixed(modelObj$fixed[[elem]])){ inits[[elem]] = modelObj$fixed[[elem]]
    }else { #use inits but replace any fixed elements with their fixed values
      inits[[elem]]=array(inits[[elem]],dim=dim(modelObj$fixed[[elem]]))
      inits[[elem]][!is.na(modelObj$fixed[[elem]])]= modelObj$fixed[[elem]][!is.na(modelObj$fixed[[elem]])]
    }
    }
  
  for(elem in c("Q","R","B","V0")) {  #if inits is a scalar to vector, make init a diagonal matrix
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
     stop(paste("Stopped in MARSSinits() because inits$",elem," is not a matrix.  See help file.\n",sep=""), call.=FALSE)
    }
    
  if( is.fixed(modelObj$fixed$x0) ){ inits$x0 = modelObj$fixed$x0 
  }else{      
    if(identical(unname(inits$x0),-99)) {  #get estimate of x0
      y1=y[,1,drop=FALSE]
      if(is.na(miss.value)){ #replace NAs with 0s
        y1[ is.na(y1) ] = 0 
      }else{ y1[ y1==miss.value ] = 0 }
      #the following is by solving for x1 using y1=Z*(D*pipi+f)+a
      pipi = solve(t(D)%*%D)%*%t(D)%*%(solve(t(inits$Z)%*%inits$Z)%*%t(inits$Z)%*%(y1-inits$A) - f)
      inits$x0 = D%*%pipi+f
    }else{ #use inits but replace any fixed elements with their fixed values
      inits[[elem]]=array(inits[[elem]],dim=dim(modelObj$fixed[[elem]]))
      inits[[elem]][!is.na(modelObj$fixed[[elem]])]= modelObj$fixed[[elem]][!is.na(modelObj$fixed[[elem]])]
    }
  }

  inits
}
