MARSSapplynames=function(obj, Y.names=NA, X.names=NA, x0.names=NA, V0.names=NA, U.names=NA, A.names=NA, R.names=NA, Q.names=NA, B.names=NA, Z.names=NA, rows=TRUE, cols=TRUE){
## Helper function to put names on the elements in an object
  theclass="marssMLE"
  if(!(class(obj) %in% c("marssMLE", "marssm"))) stop("MARSSapplynames: this function is for marssMLE and marssm objects only")
  if(class(obj)=="marssm") {obj = list(model=obj); theclass="marssm" }
  y = obj$model$data
  n = dim(y)[1]; m = dim(as.matrix(obj$model$fixed$Q))[1]
  tmp = dim(obj$model$fixed$Z)
  one.x.to.one.y = ifelse(tmp[1] == tmp[2], TRUE, FALSE)

  matnames=list()
  if(identical(Y.names,NA)) {
    if(is.null(rownames(y))) Y.names = paste("Y",seq(1, n),sep="")
    else Y.names = rownames(y) 
  }
  if(identical(X.names,NA)) X.names=paste("X",1:m,sep="")
  if(identical(x0.names,NA)) {
    tmp=obj$model$free$x0; tmp[is.na(tmp)]="fix"
    matnames$x0=paste(X.names, tmp, sep=":")
    } 
  if(identical(U.names,NA)) {
    tmp=obj$model$free$U; tmp[is.na(tmp)]="fix"
    matnames$U=paste(X.names, tmp, sep=":")
    }   
  if(identical(A.names,NA)) {
      tmp=obj$model$free$A; tmp[is.na(tmp)]="fix"
      matnames$A=paste(Y.names, tmp, sep=":")
    }         
  if(identical(R.names,NA)) {
    tmp=takediag(obj$model$free$R); tmp[is.na(tmp)]="fix"
    matnames$R=paste(Y.names, tmp, sep=":")
    }    
  if(identical(Q.names,NA)) {
    tmp=takediag(obj$model$free$Q); tmp[is.na(tmp)]="fix"
    matnames$Q=paste(X.names, tmp, sep=":")
    }      
  if(identical(B.names,NA)) {
    tmp=takediag(obj$model$free$B); tmp[is.na(tmp)]="fix"
    matnames$B=paste(X.names, tmp, sep=":")
    }
  if(identical(V0.names,NA)) matnames$V0=X.names        
  if(identical(Z.names,NA)){matnames$Z=Y.names
  }else matnames$Z=Z.names$rows
    
if(rows) {
  for(elem in c("U","A","x0","V0","Q","R","B","Z")){
    tmp=matnames[[elem]]
    if(!is.null(obj$par[[elem]])) rownames(obj$par[[elem]]) = tmp
    if(!is.null(obj$start[[elem]])) rownames(obj$start[[elem]]) = tmp
    if(!is.null(obj$model$fixed[[elem]])) rownames(obj$model$fixed[[elem]]) = tmp
    if(!is.null(obj$model$free[[elem]])) rownames(obj$model$free[[elem]]) = tmp
  }
  if(!is.null(obj$kf$xtT)) rownames(obj$kf$xtT) =  X.names
  if(!is.null(obj$states.se)) rownames(obj$states.se) =  X.names
  if(!is.null(obj$model$data) && is.null(rownames(obj$model$data))) rownames(obj$model$data) =  Y.names
}
  if(identical(Z.names,NA)){matnames$Z=X.names
  }else matnames$Z=Z.names$cols
if(cols) {
  for(elem in c("V0","Q","R","B","Z")){
    tmp=matnames[[elem]]
    if(!is.null(obj$par[[elem]])) colnames(obj$par[[elem]]) = tmp
    if(!is.null(obj$start[[elem]])) colnames(obj$start[[elem]]) = tmp
    if(!is.null(obj$model$fixed[[elem]])) colnames(obj$model$fixed[[elem]]) = tmp
    if(!is.null(obj$model$free[[elem]])) colnames(obj$model$free[[elem]]) = tmp
  }
}
  if(theclass=="marssm") return(obj$model)
  else return(obj)
  
}
