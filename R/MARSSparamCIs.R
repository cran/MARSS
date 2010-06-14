#######################################################################################################
#   MARSSparamCIs function
#   This returns CIs for ML parameter estimates
#######################################################################################################
MARSSparamCIs = function(MLEobj, method="hessian", alpha=0.05, nboot=1000) {
#this function expects a marssMLE object  as output by KalmanEM
if(!(method %in% c("hessian","parametric","innovations"))) stop("MARSSparamCIs: Current methods are hessian, parametric, and innovations")
if(!is.marssMLE(MLEobj)) stop("MARSSparamCIs needs a valid marss MLE object (this is for MLE not Bayes objects).")
if(method=="hessian")  {
    #if the model has no Hessian specified, then run emHessian to get it
    if(is.null(MLEobj$Hessian)) MLEobj = MARSShessian(MLEobj)
    #standard errors
    stderr=try(sqrt(diag(solve(MLEobj$Hessian))), silent=TRUE)      #invert the Hessian; wrap in a try() in case it fails
    if(inherits(stderr, "try-error") || any(is.na(stderr)) ) {
        paramvector = MARSSvectorizeparam(MLEobj)
        stderr=rep(NA,length(paramvector))
        warning("MARSSparamCIs: Hessian cannot be inverted; stderr = NA is being returned")
        }
    #stderr has the stderr in $par for the estimated parameters but $fixed for the fixed
    stderr.model=MARSSvectorizeparam(MLEobj, stderr)
    #need to set the se's for the fixed values at NA
    #the following loops through the elements to do this
 
    ## Order of elements is alphabetical
    free = stderr.model$model$free
    param = MLEobj$par

    parlen = 0; maxvec = NULL; par.se = NULL; par.upCI = NULL; par.lowCI = NULL
    ## Check length(parvec) matches number of free params
    for(elem in model.elem.w.V0 ) { #model.elem is speced in MARSSsettings
      free.mat = free[[elem]]
      if(any(!is.na(free.mat))) { ## if any free params
        this.mat = stderr.model$par[[elem]]
        na.mat = as.numeric(free[[elem]]) 
        na.mat[is.na(na.mat)]=NA; na.mat[!is.na(na.mat)]=0
        this.mat = this.mat + na.mat
        tmp = list(this.mat); names(tmp) = elem
	par.se = c(par.se, tmp)
	tmp = list(param[[elem]] + qnorm(1-alpha/2)*this.mat); names(tmp) = elem
        par.upCI = c(par.upCI, tmp)
	tmp = list(param[[elem]] - qnorm(1-alpha/2)*this.mat); names(tmp) = elem
        par.lowCI = c(par.lowCI, tmp)
      } #end if(any(!is.na(this.mat)))

      else { ## use free matrix since it will be all NAs (since none free)
	       par.se = c(par.se, as.numeric(free[[elem]]))
      } #else no na in free.mat
    } #end elem loop
    par.bias = NULL; par.CI.nboot = NULL
} #if method hessian
if(method %in% c("parametric", "innovations"))  {
    boot.params = MARSSboot(MLEobj, nboot=nboot, output="parameters", sim=method,
          param.gen="KalmanEM", silent=TRUE)$boot.params
    vector.lowCI = apply(boot.params,1,quantile,probs=alpha/2)
    par.lowCI = MARSSvectorizeparam(MLEobj, vector.lowCI)$par
    vector.upCI = apply(boot.params,1,quantile,probs=1-alpha/2)
    par.upCI = MARSSvectorizeparam(MLEobj, vector.upCI)$par
    vector.par.se = sqrt(apply(boot.params,1,var))
    par.se = MARSSvectorizeparam(MLEobj, vector.par.se)$par
    paramvec = MARSSvectorizeparam(MLEobj)
    vector.bias = paramvec - apply(boot.params,1,mean)
    par.bias = MARSSvectorizeparam(MLEobj, vector.bias)$par  
    par.CI.nboot = nboot
}
#set values to NA for fixed elements
for(elem in model.elem.w.V0) { #model.elem is speced in MARSSsettings
   par.se[[elem]][is.na(MLEobj$model$free[[elem]])]=NA
   par.lowCI[[elem]][is.na(MLEobj$model$free[[elem]])]=NA
   par.upCI[[elem]][is.na(MLEobj$model$free[[elem]])]=NA
   if(!is.null(par.bias)) par.bias[[elem]][is.na(MLEobj$model$free[[elem]])]=NA
   }
MLEobj$par.se = par.se
MLEobj$par.bias = par.bias
MLEobj$par.upCI = par.upCI
MLEobj$par.lowCI = par.lowCI
MLEobj$par.CI.alpha = alpha
MLEobj$par.CI.method = method
MLEobj$par.CI.nboot = par.CI.nboot
return(MLEobj)
}

