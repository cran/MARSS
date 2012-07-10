###############################################################################################################################################
#  Print method for class marssMLE. 
##############################################################################################################################################
print.marssMLE <- function (x, digits = max(3, getOption("digits")-4), ..., what="fit", silent=FALSE) 
{
if(what=="fit"){
      ## If estimated parameters present
      if (!is.null(x$par)) {
      cat("\nMARSS fit is\n")
      cat(paste("Estimation method:", x$method, "\n"))
	    if(x$control$MCInit) cat("Monte Carlo initialization with", x$controlnumInits, "random starts. \n")
      if( x$method %in% kem.methods ) 
        cat(paste("Convergence test: conv.test.slope.tol = ", x$control$conv.test.slope.tol,", abstol = ", x$control$abstol, "\n", sep=""))
	    if(x$convergence==0){
        if( x$method %in% optim.methods ){ cat("Estimation converged in", x$numIter,"iterations. \n")
	       }else{
	         if( x$numIter>x$control$minit ){ cat("Estimation converged in", x$numIter,"iterations. \n")
	         }else{ cat("Algorithm ran", x$numIter,"(=minit) iterations and convergence was reached. \n") }
	       }          
	    }
	    if ( x$convergence==1 ){
        if(x$method %in% kem.methods) tmp.msg=paste("Neither abstol nor log-log convergence test were passed.\n", sep="") else tmp.msg=""
         cat("WARNING: maxit reached at ",x$control$maxit," iter before convergence.\n", tmp.msg,
         "The likelihood and params are not at the ML values.\n",
         "Try setting control$maxit higher.\n")
      }
	    if (x$convergence==2) cat("Invalid MLE object. \n")
	    if (x$convergence==3) cat("No convergence information. \n")
      if (x$method %in% kem.methods && x$convergence %in% c(10,12)){
        if(x$convergence==10) tmp.msg=paste("maxit reached at ",x$control$maxit," iter before log-log convergence.\n", sep="") 
        if(x$convergence==12) tmp.msg=paste("maxit (=",x$control$maxit,") < min.iter.conv.test (=",x$control$min.iter.conv.test,") therefore no log-log test info.\n", sep="")           
         cat("WARNING: abstol convergence only no log-log convergence.\n", tmp.msg,
         "The likelihood and params might not be at the ML values.\n",
         "Try setting control$maxit higher.\n")
      }
      if (x$method %in% kem.methods && x$convergence==11){
         tmp.msg=paste("maxit reached at ",x$control$maxit," iter before abstol convergence.\n", sep="")
         cat("WARNING: log-log convergence only no abstol convergence.\n", tmp.msg,
         "The likelihood and params might not be at the ML values.\n",
         "Try setting control$maxit higher.\n")
      }
      if (x$method %in% kem.methods && x$convergence==13){
         tmp.msg=paste("maxit (=",x$control$maxit,") < min.iter.conv.test (=",x$control$min.iter.conv.test,") therefore no log-log test info.\n", sep="")
         cat("WARNING: maxit reached before abstol convergence and no log-log test info.\n", tmp.msg,
         "The likelihood and params are not at the ML values.\n",
         "Try setting control$maxit higher.\n")
      }
      if (x$method %in% optim.methods && x$convergence==10) 
         cat("WARNING: degeneracy of the Nelder-Mead simplex. \n")
	    if (x$method %in% kem.methods && x$convergence%in%c(52,62)) 
         cat("WARNING: Estimation was stopped at iteration", x$numIter, "due to errors.\n", 
         "Parameter estimates are those at the last iteration before stopping\n",
         "see $errors element of marssMLE object to view errors\n")
	    if (x$method %in% optim.methods && (x$convergence==51 || x$convergence==52)) 
         cat("WARNING: error or warning from the L-BFGS-B method; see component message for details\n")
      if (!(x$convergence %in% c(0,1,2,3,10,11,12,13,51,52,62)))
         cat("WARNING: convergence test errors\n")

	cat("Log-likelihood:", x$logLik, "\n")
	if(!is.null(x$AIC)) cat("AIC:", x$AIC,"  ") 
	if(!is.null(x$AICc)) cat("AICc:", x$AICc,"  ") 
	if(!is.null(x$AICbb)) cat("AICbb(innov):", x$AICbb, "  ") 
	if(!is.null(x$AICbp)) cat("AICbp(param):", x$AICbp, "  ")
  cat("\n \n") 

	cmat = as.matrix(MARSSvectorizeparam(x))
	colnames(cmat) = "Estimate"
	if(!is.null(x$par.se)) {   
	  cmat = cbind(cmat, x$par.se)
	  colnames(cmat) = c("ML.Est", "Std.Err")
	}

  if(!is.null(x$par.lowCI) & !is.null(x$par.upCI)) {  
	  cmat = cbind(cmat, x$par.lowCI)
	  colnames(cmat)[dim(cmat)[2]]="low.CI"
	  cmat = cbind(cmat, x$par.upCI)
	  colnames(cmat)[dim(cmat)[2]]="up.CI"  
	}
	if(!is.null(x$par.bias)) {  
	  cmat = cbind(cmat, x$par.bias)
	  colnames(cmat)[dim(cmat)[2]]="Est.Bias"
    cmat = cbind(cmat, as.matrix(MARSSvectorizeparam(x))+x$par.bias)
	  colnames(cmat)[dim(cmat)[2]]="Unbias.Est"
	}
	#printCoefmat(cmat, digits=digits)
	print(cmat, digits=digits)
      
  cat("\n")
  if(is.null(x$par.se)) cat("Standard errors have not been calculated. \n")
  if(!is.null(x$par.lowCI) & !is.null(x$par.upCI)){cat(paste("CIs calculated at alpha = ", x$par.CI.alpha, " via method=", x$par.CI.method, sep=""), "\n")
  }else cat("Use MARSSparamCIs to compute CIs and bias estimates.\n")
	if(!is.null(x$par.bias)){cat(paste("Bias calculated via",x$par.CI.method,"bootstrapping with",x$par.CI.nboot,"bootstraps. \n"))
  }
  
  if(!is.null(x$errors)) cat(x$errors)
   
  }else {
    cat("marssMLE object $par element is NULL.  Parameters have not been estimated.\n")
    if (x$convergence==2){ 
         cat("ERROR: marssMLE object is inconsistent or incomplete.\n No model was fit.\n marssMLE object $convergence arg = 2\n")
         }
    if (x$convergence==52 || x$convergence==62) 
         cat("WARNING: Estimation was stopped due to errors.\n", 
         "see $errors element of marssMLE object to view errors\n")
    if (x$convergence==53){ 
         cat("ERROR: Estimation was stopped in optim() due to errors returned by MARSSkf.\n", 
         "No parameter estimates are available.\n marssMLE object $convergence arg = 53\n")
         }
    }  
    cat("\n")
    return.obj=x
    }
    if(what=="model"){
       if(is.null(x$call)){ print(x$model) 
       }else{
         if(x$call$form=="base"){print(x$model)
         }else{
           print.fun=paste("print.model.",x$call$form,sep="")
           if(!exists(print.fun)){ print(x$model)
           }else eval(call(print.fun, x))
         }
        }
      return.obj=x$model
    }
    if(what=="par"){
      if(!silent) print(MARSSvectorizeparam(x))
      return.obj=MARSSvectorizeparam(x)
    }
    if(what=="xtT" | what=="states"){
      if(!silent) print(x$states)
      return.obj=x$states
    }
    if(what=="data"){
      if(!silent) print(x$model$data)
      return.obj=x$model$data
    }
    if(what=="ytT"){
      if(is.null(x$Ey)){ ytT=MARSShatyt( x )$ytT }else{ ytT=x$Ey$ytT }
      if(!silent) print(ytT)
      return.obj=ytT
    }
    if(what=="states.se"){
      if(!silent) print(x$states.se)
      return.obj=x$states.se
    }
    if(what=="states.cis"){
      if(!silent) cat("Approximate 95% confidence intervals for the states using 1.95*states.se.\n")
      rpt.list=list(up95CI=x$states+qnorm(0.975)*x$states.se,est=x$states, low95CI=x$states-qnorm(0.975)*x$states.se)
      if(!silent) print(rpt.list)
      return.obj=rpt.list
    }
    if(what %in% names(x$model$fixed)){
      the.par = parmat(x,what,t=1:max(dim(x$model$fixed[[what]])[3],dim(x$model$free[[what]])[3]))[[what]]
      if(!silent) print(the.par)
      return.obj=the.par
    }
    invisible(return.obj)  
    }
    