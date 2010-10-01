#######################################################################################################
#   MARSShessian functions
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#######################################################################################################
MARSShessian = function(MLEobj) {

  paramvector = MARSSvectorizeparam(MLEobj)

  kfNLL=function(paramvec, MLEobj){
    new.model = MARSSvectorizeparam(MLEobj, parvec=paramvec)
    y = new.model$model$data
    kf = MARSSkf(y, parList=new.model$par, missing.matrix=new.model$model$M, miss.value=new.model$model$miss.value)
    return(-kf$logLik)
  }

  #Hessian and gradient
  emhess = fdHess(paramvector, function(paramvector, MLEobj) kfNLL(paramvector, MLEobj), MLEobj)
  MLEobj$Hessian = emhess$Hessian
  MLEobj$gradient = emhess$gradient

  parSigma = try(solve(MLEobj$Hessian), silent=TRUE)
  if(inherits(parSigma, "try-error")) {
    warning("MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix")
    parSigma=NULL
  }
  MLEobj$parSigma = parSigma
  MLEobj$parMean = paramvector

  return(MLEobj)
}
