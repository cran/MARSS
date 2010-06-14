###############################################################################################################################################
#  Print method for class marssMLE. 
###############################################################################################################################################

print.marssMLE <- function (x, digits = max(3, getOption("digits")-4), ...) 
    {
      ## If estimated parameters present
      if (!is.null(x$par)) {
        cat("\nMARSS fit is\n")
        cat(paste("Estimation method: ", x$method, "\n"))
	if(x$control$MCInit) cat("Monte Carlo initialization with", x$controlnumInits, "random starts. \n")

      if (x$convergence==1) cat("WARNING: Algorithm reached",x$control$maxit,"iterations before hitting tolerance threshold. \n")
	    if (x$convergence==0) cat("Estimation converged in", x$numIter,"iterations. \n")
	    if (x$convergence==52) 
         cat("WARNING: Estimation was stopped at iteration", x$numIter, "due to errors.\n", 
         "Parameter estimates are those at the last iteration before stopping\n",
         "see $errors element of marssMLE object to view errors\n")
	    if (x$convergence==10) 
         cat("WARNING: Tolerance reached at iter=",x$numIter," before param ests converged.\n",  
         "The likelihood and params are not at the ML values.\n",
         "Try setting control$minit higher or control$abstol lower\n")

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
    if (x$convergence==52) 
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
  
