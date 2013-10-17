#######################################################################################################
#   MARSShessian functions
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#######################################################################################################
MARSShessian = function(MLEobj) {

  fun="MARSSkf"

  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MLEobj[["marss"]][["data"]] #must have time going across columns
  marss.object = MLEobj[["marss"]]
  free = marss.object[["free"]]
  fixed = marss.object[["fixed"]]
  par.dims=attr(marss.object,"model.dims")
  tmp.par=MLEobj[["par"]]
  
  #The code is used to set things up to use MARSSvectorizeparam to just select inits for the estimated parameters
  tmp.MLEobj = MLEobj
  for(elem in c("Q","R","V0")){ #need the chol for these
    #I think you can have time-varying but I need to figure how to compute fixed and free using the t=1 par; 
    #and I probably have to put constraints on how it is time-varying
    if(par.dims[[elem]][3] != 1) stop("Stopped in MARSShessian(): Time-varying variances are not allowed in the Hessian computation.")
    f=sub3D(fixed[[elem]],t=1) #not time-varying by constraint above
    d=sub3D(free[[elem]],t=1) #free[[elem]] is required to be time constant
    the.par=orig.par=unvec(f+d%*%tmp.par[[elem]], dim=par.dims[[elem]][1:2])
    is.zero=diag(orig.par)==0   #where the 0s on diagonal are
    if(any(is.zero)) diag(the.par)[is.zero]=1    #so the chol doesn't fail if there are zeros on the diagonal
    the.par=t(chol(the.par))  #transpose of chol
    if(any(is.zero)) diag(the.par)[is.zero]=0  #set back to 0
    if(!is.fixed(free[[elem]])){
      tmp.MLEobj[["par"]][[elem]] = solve(crossprod(d))%*%t(d)%*%(vec(the.par)-f) #from f+Dm=M and if f!=0, D==0
    }else{ tmp.MLEobj[["par"]][[elem]] = matrix(0,0,1) }
    #when being passed to optim, pars for var-cov mat is the chol, so need to reset free and fixed
    #step 1, compute the D matrix corresponding to upper.tri=0 at in t(chol)
    tmp.list.mat=fixed.free.to.formula(sub3D(tmp.MLEobj[["marss"]][["fixed"]][[elem]],t=1),sub3D(tmp.MLEobj[["marss"]][["free"]][[elem]],t=1),par.dims[[elem]][1:2])
    tmp.list.mat[upper.tri(tmp.list.mat)]=0   #set upper tri to zero
    tmp.MLEobj[["marss"]][["free"]][[elem]]=convert.model.mat(tmp.list.mat)[["free"]]
    #step 2, set the fixed part to the t(chol); need to take the chol of fixed elements; upper tri will be set to 0
    tmp.fixed=unvec(f,dim=par.dims[[elem]][1:2])  #by definition the estimated elements will have f=0 since this is a varcov mat
    is.zero = diag(tmp.fixed)==0 #need to deal with zeros on diagonal
    if(any(is.zero)) diag(tmp.fixed)[is.zero]=1    #so the chol doesn't fail if there are zeros on the diagonal
    chol.fixed=t(chol(tmp.fixed)) #chol of the fixed part of var-cov matrix
    if(any(is.zero)) diag(chol.fixed)[is.zero]=0  #reset back to zero   
    tmp.MLEobj[["marss"]][["fixed"]][[elem]][,,1]=vec(chol.fixed) #above required that dim3 of fixed$Q and R is 1
    }
  # will return the inits only for the estimated parameters
  paramvector = MARSSvectorizeparam(tmp.MLEobj)
  
  MLEobj=tmp.MLEobj
    
#   kfNLL=function(paramvec, MLEobj){
#     #neglogLik is defined in MARSSoptim
#     return(neglogLik(paramvec, MLEobj=MLEobj))
#   }
  #Hessian and gradient
  emhess = fdHess(paramvector, function(paramvector, MLEobj) kfNLL(paramvector, MLEobj), MLEobj)
  MLEobj$Hessian = emhess$Hessian
  rownames(MLEobj$Hessian)=names(paramvector)
  colnames(MLEobj$Hessian)=names(paramvector)
  MLEobj$gradient = emhess$gradient

  parSigma = try(solve(MLEobj$Hessian), silent=TRUE)
  if(inherits(parSigma, "try-error")) {
    warning("MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix")
    parSigma=NULL
  }
  MLEobj$parSigma = parSigma
  MLEobj$parMean = paramvector

  #This is the TRANSFORMED MLEobj
  return(MLEobj)
}

kfNLL = function(x, MLEobj=NULL){  #NULL assignment needed for optim call syntax
  #MLEobj is tmp.MLEobj so has altered free and fixed
  #x is the paramvector
  MLEobj = MARSSvectorizeparam(MLEobj, x)
  marss.object = MLEobj[["marss"]]
  free=marss.object[["free"]]
  fixed=marss.object[["fixed"]]
  pars=MLEobj[["par"]]
  par.dims=attr(marss.object,"model.dims")
  for(elem in c("Q","R","V0")){
    if(!is.fixed(free[[elem]])) #recompute par if needed since par in parlist is transformed
    {        
      d=sub3D(free[[elem]],t=1) #this will be the one with the upper tri zero-ed out but ok since symmetric
      par.dim=par.dims[[elem]][1:2]
      #t=1 since D not allowed to be time-varying?
      L=unvec(free[[elem]][,,1]%*%pars[[elem]],dim=par.dim) #this by def will have 0 row/col at the fixed values
      the.par = tcrossprod(L)#L%*%t(L)
      #from f+Dm=M and if f!=0, D==0 so can leave off f
      MLEobj[["par"]][[elem]]=solve(crossprod(d))%*%t(d)%*%vec(the.par)
      #solve(t(d)%*%d)%*%t(d)%*%vec(the.par)
    }
  } #end for over elem
  #kfsel selects the Kalman filter / smoother function based on MLEobj$fun.kf
  negLL = MARSSkf( MLEobj, only.logLik=TRUE, return.lag.one=FALSE )$logLik
  
  -1*negLL
}

MARSShessian.backtrans = function(MLEobj.hessian, par.hessian){
  #MLEobj is you original untransformed MLEobj
  #MLEobj.hessian is a transformed version with the chol transformation for variances
  #par.hessian is a vector of parameters where the variances are in chol form
  #Goal is to put the par variances values in 
  
  #first put the parameters into the MLEobj.hession object
  MLEobj.hessian=MARSSvectorizeparam(MLEobj.hessian, par.hessian)

  #chol back transformation
  par.dims=attr(MLEobj.hessian[["marss"]],"model.dims")
  for(elem in c("Q","R","V0")){   #this works because by def fixed and free blocks of var-cov mats are independent
  if(!is.fixed(MLEobj.hessian$marss$free[[elem]])) #if not estimated then there won't be a par element
  {
    d=sub3D(MLEobj.hessian[["marss"]][["free"]][[elem]],t=1) #this will be the one with the upper tri zero-ed out but ok since symmetric
    par.dim=par.dims[[elem]][1:2]
    L=unvec(MLEobj.hessian[["marss"]][["free"]][[elem]][,,1]%*%MLEobj.hessian$par[[elem]],dim=par.dim) #this by def will have 0 row/col at the fixed values
    the.par = tcrossprod(L) #L%*%t(L)
    MLEobj.hessian[["par"]][[elem]]=solve(crossprod(d))%*%t(d)%*%vec(the.par)
  }
} #end for

#now the MLEobj par elements are back transformed
return( MARSSvectorizeparam(MLEobj.hessian) )
}
