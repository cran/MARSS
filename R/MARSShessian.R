#######################################################################################################
#   MARSShessian functions
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#######################################################################################################
MARSShessian = function(MLEobj, fun="MARSSkf") {

  paramvector = MARSSvectorizeparam(MLEobj)

  kfNLL=function(paramvec, MLEobj, fun){
    new.MLEobj = MARSSvectorizeparam( MLEobj, parvec=paramvec )
    kf=eval(call(fun, new.MLEobj))
    return(-kf$logLik)
  }

  #Hessian and gradient
  emhess = fdHess(paramvector, function(paramvector, MLEobj, fun) kfNLL(paramvector, MLEobj, fun), MLEobj, fun)
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

  return(MLEobj)
}
