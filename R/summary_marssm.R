###############################################################################################################################################
#  Summary method for class marssm. 
###############################################################################################################################################

summary.marssm <- function (object, ...) 
{

   n = dim(object$data)[1]; m = dim(object$fixed$x0)[1]
   cat(paste("Model Structure is\n","m: ",m," state process(es)\n","n: ",n," observation time series\n",sep=""))

      rpt.list = list()
      en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")
      dim.tmp = list(Z=c(n,m), A=c(n,1), R=c(n,n), B=c(m,m), U=c(m,1), Q=c(m,m), x0=c(m,1), V0=c(m,m))

	xnames = paste("X",1:m,sep="")
	if( !is.null(object[["X.names"]]) ) xnames = object$X.names
	ynames = paste("Y",1:n,sep="")
	if( !is.null(rownames(object[["data"]])) ) ynames = rownames(object$data)
      
  for (elem in en) {
    #list matrix version of the model    
    Tmax=max(dim(object$fixed[[elem]])[3],dim(object$free[[elem]])[3])
    for(t in 1:Tmax){
      tmp = fixed.free.to.formula( sub3D(object$fixed[[elem]],t=t),sub3D(object$free[[elem]],t=t),dim.tmp[[elem]] )
      if(Tmax==1) rpt.list[[elem]] = tmp
      if(t==1 & Tmax>1) rpt.list[[elem]]=array(list(),dim=c(dim(tmp),Tmax))
      if(Tmax>1) rpt.list[[elem]][,,t] = tmp
    } # for t in Tmax
	  if(elem %in% c("B","Q","Z","V0")) colnames(rpt.list[[elem]])=xnames
	  if(elem %in% c("x0","U")) colnames(rpt.list[[elem]])="X"
	  if(elem %in% c("R")) colnames(rpt.list[[elem]])=ynames
	  if(elem %in% c("A")) colnames(rpt.list[[elem]])="Y"
	  if(elem %in% c("B","U","Q","x0","V0")) rownames(rpt.list[[elem]])=xnames
	  if(elem %in% c("Z","A","R")) rownames(rpt.list[[elem]])=ynames
  } #for elem
  
  print(rpt.list, quote=FALSE)
  invisible(object)
          
  }