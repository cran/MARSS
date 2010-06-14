###############################################################################################################################################
#  Summary method for class marssm. 
###############################################################################################################################################

summary.marssm <- function (object, ...) 
    {
      n = dim(object$fixed$Z)[1]; m = dim(object$fixed$Z)[2]
      cat(paste("Model Structure is\n","m: ",m," state process(es)\n","n: ",n," observation time series\n",sep=""))

      rpt.list = NULL
      en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")

	# bug 600
	gnames = 1:dim(object$fixed$Z)[2]
	if( !is.null(colnames(object$fixed$Z)) ) gnames = colnames(object$fixed$Z)
      
      for (elem in en) {

	## Mark free params with element.group
	#tmp.free =  paste("g", object$free[[elem]], sep="")
	dd = dim(object$free[[elem]])
	named.mat = matrix(gnames[ object$free[[elem]] ], nrow=dd[1], ncol=dd[2])
	tmp.free =  paste(elem, named.mat, sep=".")

	## Mark fixed params with ()
	tmp = paste("(", object$fixed[[elem]], sep="")
	tmp = paste(tmp, ")", sep="")

	## Combine free & fixed
	tmp[tmp == "(NA)"] = tmp.free[tmp == "(NA)"] 

	rpt.list[[elem]] = matrix(tmp, nrow=nrow(object$fixed[[elem]]))
      }

      print(rpt.list, quote=FALSE)
      invisible(object)

    }            
