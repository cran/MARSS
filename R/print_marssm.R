###############################################################################################################################################
#  Print MARSS model structure 
###############################################################################################################################################

print.marssm <- function (x, ...) 
    {    
      if (is.null(dim(x$data))) n = 1
      else n = dim(x$data)[1]
      m=dim(x$fixed$Q)[1]
      
      cat(paste("\nModel Structure is\n","m: ", m," state process(es)\n", "n: ", n," observation time series\n",sep=""))
      tmp = NULL

      ## Print constraints for each parameter
      rpt = describe.marssm(x)
      if(is.design(x$fixed$Z)){
        ## Print sites by group
        group.names=rpt$Z
        if(!is.null(rownames(x$data))) data.names=rownames(x$data)
        else data.names = paste("Y",1:n,sep="")
        for (i in unique(group.names)) {
	         cat(paste("State", i, ": "))
	         cat(data.names[group.names==i], "\n")
            }
        cat("\n")
        }
	  
      for (el in names(x$fixed)) {
	# if constraint is in English, print it
  cat(el, ": ", rpt[[el]], "\n")
      }  
}

