find.degenerate=function(MLEobj){
#A simple helper function to find the degenerate variance parameters
##error checking
  tmp = is.marssMLE(MLEobj)
  if(!isTRUE(tmp)) {
      stop("find.degenerate: marssMLE object is incomplete or inconsistent.\n", call.=FALSE)
    }
  if(is.null(MLEobj$par) || is.null(MLEobj$iter.record)) stop("find.degenerate: MLE object does not have par or iter.record element.\n", call.=FALSE) 
  if(MLEobj$method != "kem") stop("find.degenerate: This function is only for marssMLE objects generated from method=kem.\n", call.=FALSE) 

#find the variance components and plot them
names.iter=colnames(MLEobj$iter.record$par)
names.sub=strsplit(names.iter,"\\.")
num.names = length(names.sub)
b=NULL
for(j in 1:num.names)b=c(b,names.sub[[j]][1])
b.varcov = b[b %in% c("R","Q")]
num.varcov = length(b.varcov)
if(num.varcov<4){ par(mfrow=c(1,num.varcov))
}else par(mfrow=c(ceiling(sqrt(num.varcov)),ceiling(num.varcov/ceiling(sqrt(num.varcov)))))
for( j in 1:num.names ){
   if(b[j] %in% c("R","Q")) {
      ylims=c(min(log(abs(MLEobj$iter.record$par[,j])))-.05,max(log(abs(MLEobj$iter.record$par[,j])))+.05)
      len=length(MLEobj$iter.record$par[,j])
      test.len2=dim(MLEobj$iter.record$par)[1]
    	test.len1=max(1,test.len2-9)
      test.len=(MLEobj$numIter-min(test.len2-1,9)):MLEobj$numIter      
	test=lm(log(abs(MLEobj$iter.record$par[max(1,len-9):len,j]))~log(test.len))
	if(abs(test$coef[2])>0.5){thecol="red"}else thecol="black"
      plot(log(test.len),log(abs(MLEobj$iter.record$par[1:test.len2,j])),ylim=ylims, type="l",xlab="log(iter)",ylab="abs(log(param))",col=thecol)
	abline(v=log(test.len[1]),col="green")
	abline(v=log(test.len[10]),col="green")
	    if(MLEobj$numIter>50){ warn.txt = "degenerate"
	    }else warn.txt="not converged"
      if(abs(test$coef[2])>0.5){title(paste(names.iter[j],format(test$coef[2],digits=2),"\n",warn.txt))
      }else title(paste(names.iter[j],format(test$coef[2],digits=2)))
      }
}
}

      