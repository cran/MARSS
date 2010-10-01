###############################################################################################################################################
#  Summary method for class marssm. 
###############################################################################################################################################

summary.marssm <- function (object, ...) 
    {
      n = dim(object$fixed$Z)[1]; m = dim(object$fixed$Z)[2]
      cat(paste("Model Structure is\n","m: ",m," state process(es)\n","n: ",n," observation time series\n",sep=""))

      rpt.list = list()
      en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")

	xnames = paste("X",1:dim(object$fixed$Z)[2],sep="")
	if( !is.null(colnames(object$fixed$Z)) ) xnames = colnames(object$fixed$Z)
	ynames = paste("Y",1:dim(object$fixed$Z)[1],sep="")
	if( !is.null(rownames(object$fixed$Z)) ) ynames = rownames(object$fixed$Z)
      
      for (elem in en) {

	## Mark free params with element.group
	tmp.free =  paste(elem, object$free[[elem]], sep=".")

	## Mark fixed params with ()
	tmp = paste("(", object$fixed[[elem]], sep="")
	tmp = paste(tmp, ")", sep="")

	## Combine free & fixed
	tmp[tmp == "(NA)"] = tmp.free[tmp == "(NA)"] 

	rpt.list[[elem]] = matrix(tmp, nrow=nrow(object$fixed[[elem]]))
	if(elem %in% c("B","Q","Z","V0")) colnames(rpt.list[[elem]])=xnames
	if(elem %in% c("x0","U")) colnames(rpt.list[[elem]])="X"
	if(elem %in% c("R")) colnames(rpt.list[[elem]])=ynames
	if(elem %in% c("A")) colnames(rpt.list[[elem]])="Y"
	if(elem %in% c("B","U","Q","x0","V0")) rownames(rpt.list[[elem]])=xnames
	if(elem %in% c("Z","A","R")) rownames(rpt.list[[elem]])=ynames
  }
  
  print(rpt.list, quote=FALSE)
  invisible(object)
          
  }