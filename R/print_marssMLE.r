###############################################################################################################################################
#  Print method for class marssMLE. 
##############################################################################################################################################
print.marssMLE <- function (x, digits = max(3, getOption("digits")-4), ...) 
    {
      ## If estimated parameters present
      if (!is.null(x$par)) {
      cat("\nMARSS fit is\n")
      cat(paste("Estimation method:", x$method, "\n"))
	    if(x$control$MCInit) cat("Monte Carlo initialization with", x$controlnumInits, "random starts. \n")
      if(is.null(x$control$abstol)){ cat(paste("Convergence test: log-log plot with slope tol of", x$control$conv.test.slope.tol, "\n"))
      }else cat(paste("Convergence test: delta logLik with abstol of", x$control$abstol, "\n"))
	    if (x$convergence==0){
        if( x$method %in% optim.methods || x$numIter>x$control$minit ){ cat("Estimation converged in", x$numIter,"iterations. \n")
	       }else cat("Algorithm ran", x$numIter,"(=minit) iterations and convergence was reached. \n")
	    }
      if (x$convergence==1) cat("WARNING: Algorithm reached",x$control$maxit,"iterations before hitting abstol threshold. \n")
	    if (x$convergence==2) cat("Invalid MLE object. \n")
	    if (x$convergence==3) cat("MLE object computed via", x$numIter, "iterations but no convergence information. \n")
      if (x$method %in% kem.methods && x$convergence==10) 
         cat("WARNING: maxit reached,",x$control$maxit,"iter, before convergence. \n",
         "The likelihood and params are not at the ML values.\n",
         "Try setting control$maxit higher.\n")
      if (x$method %in% optim.methods && x$convergence==10) 
         cat("WARNING: degeneracy of the Nelder-Mead simplex. \n")
      if (x$convergence==11){ 
         if( x$numIter>x$control$minit ){ 
            cat("WARNING: abstol reached at iter",x$numIter,"before convergence(log-log test).\n",  
            "The likelihood and params are not at the ML values.\n",
            "Try setting control$minit higher.\n")
	       }else 
            cat("WARNING: Algorithm ran for",x$numIter,"(=minit) iterations but did not converge(log-log test). \n",
            "The likelihood and params are not at the ML values.\n",
            "Try setting control$minit higher or don't use abstol.\n",
            "fyi: abstol was reached and algorithm stopped when minit reached.\n")
         }
	    if (x$method %in% kem.methods && (x$convergence==52 || x$convergence==62)) 
         cat("WARNING: Estimation was stopped at iteration", x$numIter, "due to errors.\n", 
         "Parameter estimates are those at the last iteration before stopping\n",
         "see $errors element of marssMLE object to view errors\n")
	    if (x$method %in% optim.methods && (x$convergence==51 || x$convergence==52)) 
         cat("WARNING: error or warning from the L-BFGS-B method; see component message for details\n")
      if (!(x$convergence %in% c(0,1,2,3,10,11,51,52,62)))
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
	  x.tmp = x; x.tmp$par = x$par.se 
	  cmat = cbind(cmat, MARSSvectorizeparam(x.tmp))
	  colnames(cmat) = c("ML.Est", "Std.Err")
	}

  if(!is.null(x$par.lowCI) & !is.null(x$par.upCI)) {  
	  x.tmp = x; x.tmp$par = x$par.lowCI 
	  cmat = cbind(cmat, MARSSvectorizeparam(x.tmp))
	  colnames(cmat)[dim(cmat)[2]]="low.CI"
	  x.tmp = x; x.tmp$par = x$par.upCI 
	  cmat = cbind(cmat, MARSSvectorizeparam(x.tmp))
	  colnames(cmat)[dim(cmat)[2]]="up.CI"  
	}
	if(!is.null(x$par.bias)) {  
	  x.tmp = x; x.tmp$par = x$par.bias 
	  cmat = cbind(cmat, MARSSvectorizeparam(x.tmp))
	  colnames(cmat)[dim(cmat)[2]]="Est.Bias"
    cmat = cbind(cmat, as.matrix(MARSSvectorizeparam(x))+MARSSvectorizeparam(x.tmp))
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
    invisible(x) 
    }