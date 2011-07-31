MARSSLLprofile = function(MLEobj, param=NULL, x=NULL, LLlim=3, pstep=0.01, max.steps=20, plot=TRUE){
default.step=0.01
if(missing(MLEobj) || !is.marssMLE(MLEobj))
  stop("Stopped in MARSSLLprofile: you need to pass in a valid MLEobj",.call=FALSE)
if(is.null(param)){
  param=names(MARSSvectorizeparam(MLEobj))
  }
if(is.null(x) && is.null(LLlim) )
  stop("Stopped in MARSSLLprofile: either x or LLlim must be passed in.",.call=FALSE)
if(!is.null(LLlim) && (!is.numeric(LLlim) || length(LLlim)!=1) )
  stop("Stopped in MARSSLLprofile: LLlim needs to be numeric and length 1.",.call=FALSE)
if(!is.null(LLlim) ) LLlim=abs(LLlim)
if(!is.null(x) && (!is.vector(x) || length(x)<1 || !is.numeric(x) ) )
  stop("Stopped in MARSSLLprofile: x must be a vector of numeric values",.call=FALSE)
if(!all(param %in% names(MARSSvectorizeparam(MLEobj))))
  stop("Stopped in MARSSLLprofile: One of the param names does not match the free param names.\nType the MLEobj at the command line to see parameter names.",.call=FALSE)
method=MLEobj$method
tmp=list(); pnames=names(MLEobj$par)
if(!is.list(pstep)){tmp[pnames]=pstep; pstep=tmp
}else{ tmp=pstep; tmp[pnames[!(pnames %in% names(pstep))]]=default.step; pstep=tmp }
tmp=list()
if(!is.list(x)){tmp[pnames]=list(x); x=tmp
  }else{ tmp=x; tmp[pnames[!(pnames %in% names(x))]]=NULL; x=tmp }
if(plot){
  nr=ceiling(sqrt(length(param)))
  nc=ceiling(length(param)/nr)
  par(mfrow=c(nr,nc))
  }
rtn.list=list()
for(pname in param){
  tmp.MLEobj=MLEobj
  tmp.MLEobj$control$silent=TRUE
  tmp.MLEobj$method=method
  profLL=c()
  el=unlist(strsplit(pname,".",fixed=TRUE))[1]  #the bit before the . is the parameter name
  if(!(el %in% names(MLEobj$par))) next
  gname = substr(pname,nchar(el)+2,nchar(pname))
  if(!any(tmp.MLEobj$model$free[[el]]==gname)) next
  tmp.MLEobj$model$fixed[[el]][which(MLEobj$model$free[[el]]==gname)]=MLEobj$par[[el]][which(MLEobj$model$free[[el]]==gname)]
  tmp.MLEobj$model$free[[el]][which(MLEobj$model$free[[el]]==gname)]=NA
  if(!is.marssm(tmp.MLEobj$model)) next
  if(!is.null(x[[el]])){
     prange=x[[el]]
     for(p in prange){
      tmp.MLEobj$model$fixed[[el]][which(MLEobj$model$free[[el]]==gname)]=p
      if(method=="kem") this.LL=MARSSkem(tmp.MLEobj)$logLik 
      if(method=="BFGS") this.LL=MARSSoptim(tmp.MLEobj)$logLik 
      profLL=c(profLL,this.LL)
     }
  }else { #use LLlim
   j=1; prange=c()
   p=MLEobj$par[[el]][which(MLEobj$model$free[[el]]==gname)][1]
   while(j<=max.steps){
     prange=c(prange,p)
     tmp.MLEobj$model$fixed[[el]][which(MLEobj$model$free[[el]]==gname)]=p
     if(method=="kem") this.LL=MARSSkem(tmp.MLEobj)$logLik 
     if(method=="BFGS") this.LL=MARSSoptim(tmp.MLEobj)$logLik 
     profLL=c(profLL,this.LL)
     if(this.LL<(MLEobj$logLik-abs(LLlim))) break
     if(el %in% c("Q","R","V0")){p=p*exp(-abs(pstep[[el]]))}else p=p-abs(pstep[[el]])
     j=j+1   
   }
   j=1
   p=MLEobj$par[[el]][which(MLEobj$model$free[[el]]==gname)][1]
   while(j<=max.steps){
     if(el %in% c("Q","R","V0")){p=p*exp(abs(pstep[[el]]))}else p=p+abs(pstep[[el]])
     prange=c(prange,p)
     tmp.MLEobj$model$fixed[[el]][which(MLEobj$model$free[[el]]==gname)]=p
     if(method=="kem") this.LL=MARSSkem(tmp.MLEobj)$logLik 
     if(method=="BFGS") this.LL=MARSSoptim(tmp.MLEobj)$logLik 
     profLL=c(profLL,this.LL)
     if(this.LL<(MLEobj$logLik-abs(LLlim))) break
     j=j+1   
   }
  }
  if(plot){
    if(el %in% c("Q","R","V0")){ x.plot=log(prange) }else x.plot=prange
    plot(sort(x.plot),profLL[sort(x.plot,index.return=TRUE)$ix],type="l",ylab="LL",xlab="Value")
    v=MLEobj$par[[el]][which(MLEobj$model$free[[el]]==gname)][1]
    if(el %in% c("Q","R","V0")) v=log(v) 
    abline(v=v)
    abline(h=max(profLL)-1.92,col="red")
    title(main=pname)
    }
  rtn.list[[pname]]=cbind(prange,profLL) 
}
return(rtn.list)
}