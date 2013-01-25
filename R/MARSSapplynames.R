MARSSapplynames=function(obj, Y.names=NA, X.names=NA){
## Helper function to put names on the elements in an object
  theclass="marssMLE"
  if(!(class(obj) %in% c("marssMLE", "marssm")))
     stop("Stopped in MARSSapplynames() because this function is for marssMLE and marssm objects only.\n", call.=FALSE)
  if(class(obj)=="marssm") {obj = list(model=obj); theclass="marssm" }
  y = obj$model$data
  n = dim(y)[1]; m = dim(obj$model$fixed$x0)[1]
  one.x.to.one.y = ifelse(n == m, TRUE, FALSE)

  if(identical(Y.names,NA)) {
    if(is.null(rownames(y))) Y.names = paste("Y",seq(1, n),sep="")
    else Y.names = rownames(y) 
  }
  if(identical(X.names,NA)) {
    if(is.null(obj[["model"]][["X.names"]])) X.names = paste("X",seq(1, m),sep="")
    else X.names = obj$model$X.names 
  }
  if(length(Y.names)!=n) stop("Stopped in MARSSapplynames() because length of Y.names is not n.\n", call.=FALSE)
  if(length(X.names)!=m) stop("Stopped in MARSSapplynames() because length of X.names is not m.\n", call.=FALSE)
  
  matnames=list()
  model.elem=names(obj$model$fixed)
  for(elem in model.elem){
    if(is.null(colnames(obj[["model"]][["free"]][[elem]])) && dim(obj$model$free[[elem]])[2]>0){ 
      colnames(obj$model$free[[elem]]) = paste(seq(1, dim(obj$model$free[[elem]])[2]),sep="") #paste(elem,seq(1, dim(obj$model$free[[elem]])[2]),sep="")
    }    
    matnames[[elem]] = colnames(obj$model$free[[elem]])
  }       
    
  for(elem in model.elem){
    if(!is.null(obj[["par"]][[elem]]) & is.null(rownames(obj$par[[elem]]))) rownames(obj$par[[elem]]) = colnames(obj$model$free[[elem]])
  }
  if(!is.null(obj[["kf"]][["xtT"]])) rownames(obj$kf$xtT) =  X.names
  if(!is.null(obj[["states.se"]])) rownames(obj$states.se) =  X.names
  if(!is.null(obj[["model"]][["data"]]) && is.null(rownames(obj$model$data))) rownames(obj$model$data) =  Y.names

  if(theclass=="marssm") return(obj$model)
  else return(obj)
  
}
