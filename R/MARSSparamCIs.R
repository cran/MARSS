#######################################################################################################
#   MARSSparamCIs function
#   This returns CIs for ML parameter estimates
#######################################################################################################
MARSSparamCIs = function(MLEobj, method="hessian", alpha=0.05, nboot=1000) {
#this function expects a marssMLE object
#it will add some vectors of standard errors, biases, low/up CIs to MLEobj
if(!(method %in% c("hessian","parametric","innovations"))) stop("Stopped in MARSSparamCIs(). Current methods are hessian, parametric, and innovations.\n", call.=FALSE)
if(!is.marssMLE(MLEobj)) stop("Stopped in MARSSparamCIs(). This function needs a valid marss MLE object.\n", call.=FALSE)
paramvec = MARSSvectorizeparam(MLEobj)
if(length(paramvec)==0) stop("Stopped in MARSSparamCIs(). No estimated parameter elements.\n", call.=FALSE)
paramnames=names(paramvec)
if(method=="hessian")  {
    #if the model has no Hessian specified, then run emHessian to get it
    if(is.null(MLEobj[["Hessian"]])) MLEobj.hessian = MARSShessian(MLEobj)
    #standard errors
    hessInv = try(solve(MLEobj.hessian$Hessian))
    vector.par.se =rep(NA,length(paramvec))
    if(inherits(hessInv, "try-error") ) {
        warning("MARSSparamCIs: Hessian cannot be inverted; stderr = NA is being returned")
    }else{
    is.neg=diag(hessInv)<0
    vector.par.se[!is.neg] = try(sqrt(diag(hessInv)[!is.neg]), silent=TRUE)      #invert the Hessian; wrap in a try() in case it fails
    if(inherits(vector.par.se, "try-error") || any(is.nan(vector.par.se))) {
        warning("MARSSparamCIs: Some of the hessian diagonals cannot be square-rooted; stderr = NA is being returned")
        }
    }
    paramvec.hessian =  MARSSvectorizeparam(MLEobj.hessian)
    #this is in Cholessky transformed form
    vector.par.upCI = paramvec.hessian + qnorm(1-alpha/2)*vector.par.se
    vector.par.lowCI = paramvec.hessian - qnorm(1-alpha/2)*vector.par.se
    #need to back transform
    vector.par.upCI = MARSShessian.backtrans(MLEobj.hessian, vector.par.upCI)
    vector.par.lowCI = MARSShessian.backtrans(MLEobj.hessian, vector.par.lowCI)
    vector.par.bias = NULL; par.CI.nboot = NULL
} #if method hessian

if(method %in% c("parametric", "innovations"))  {
    boot.params = MARSSboot(MLEobj, nboot=nboot, output="parameters", sim=method,
          param.gen="MLE", silent=TRUE)$boot.params
    vector.par.lowCI = apply(boot.params,1,quantile,probs=alpha/2)
    vector.par.upCI = apply(boot.params,1,quantile,probs=1-alpha/2)
    vector.par.se = sqrt(apply(boot.params,1,var))
    paramvec = MARSSvectorizeparam(MLEobj)
    vector.par.bias = paramvec - apply(boot.params,1,mean)
    par.CI.nboot = nboot
}
names(vector.par.se)=paramnames
if(!is.null(vector.par.bias)) names(vector.par.bias)=paramnames
names(vector.par.upCI)=paramnames
names(vector.par.lowCI)=paramnames
for(val in c("par.se","par.bias","par.upCI","par.lowCI")){
  tmp.vec=get(paste("vector.",val,sep=""))
  if(!is.null(tmp.vec)) MLEobj[[val]] = MARSSvectorizeparam(MLEobj,tmp.vec)[[val]]
  else MLEobj[[val]] = NULL
}
MLEobj$par.se = MARSSvectorizeparam(MLEobj,vector.par.se)[["par"]]
if(!is.null(vector.par.bias)) MLEobj$par.bias = MARSSvectorizeparam(MLEobj,vector.par.bias)[["par"]]
  else MLEobj$par.bias = NULL
MLEobj$par.upCI = MARSSvectorizeparam(MLEobj,vector.par.upCI)[["par"]]
MLEobj$par.lowCI = MARSSvectorizeparam(MLEobj,vector.par.lowCI)[["par"]]
MLEobj$par.CI.alpha = alpha
MLEobj$par.CI.method = method
MLEobj$par.CI.nboot = par.CI.nboot
return(MLEobj)
}

