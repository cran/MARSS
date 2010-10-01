#####
#
# is.marssMLE()
# Check that the marssMLE object has all the parts it needs 
# and that these have the proper size and form
#
#####

is.marssMLE <- function(MLEobj) 
{
  if(class(MLEobj) != "marssMLE") stop("Stopped in is.marssMLE() because object class is not marssMLE.\n", call.=FALSE)
  
  msg = c()
 ## Check for required components
  el = c("model", "start", "control", "method")
  if( !all(el %in% names(MLEobj)) ){
    msg = c(msg, paste("Element", el[!(el %in% names(MLEobj))], "is missing from object.\n"))

  ## If required components present, check that the model is valid
  }else msg = is.marssm(MLEobj$model)        #returns TRUE or a vector of msgs

  ## is.wholenumber() borrowed from is.integer example
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

  ## If model is OK, check inits and params
  if (isTRUE(msg)) {

  msg = c()

  ## Set m, n, and data
  n = dim(MLEobj$model$fixed$Z)[1]
  m = dim(MLEobj$model$fixed$Z)[2]
  dat = MLEobj$model$data
  
  ## bug 584: check for T=1
  if( !is.null(dim(dat)) && !is.numeric(dat) ) msg = c(msg, "Data must be numeric.\n")
  if( !is.null(dim(dat)) && dim(dat)[2] == 1 ) msg = c(msg, "Data has only one time point.\n")

  ## Set up element names
  en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")
  
  ## Correct dimensions for reporting
  correct.dim1 = c(n,n,n,m,m,m,m,m)
  correct.dim2 = c(m,1,n,m,1,m,1,m)
  names(correct.dim1) = names(correct.dim2) = en

  ## Check initial values consistency 
  init.null = dim.init = NULL

  for (el in en) {
    init.null.flag <- ( is.null(MLEobj$start[[el]]) || !is.numeric(MLEobj$start[[el]]) )
    dim.init.flag = FALSE
 
    if (!init.null.flag) { 
      dim.init.flag <- MARSScheckdims(el, MLEobj$start, n, m) 
    }
 
    init.null <- c(init.null, init.null.flag)
    dim.init <- c(dim.init, dim.init.flag)  
  }
    problem <- any(c(init.null, dim.init))

  if (problem) {    
    if(any(init.null)) {
      msg = c(msg, paste("Missing or non-numeric initial value", en[init.null],"\n"))
    }
    if(any(dim.init)) {
      msg = c(msg, paste("Dims for initial value ", en[dim.init], " do not match Z. Dims should be ", correct.dim1[dim.init], " x ", correct.dim2[dim.init],"based on Z dims.\n", sep=""))
    }    
  }

  ## Check params consistency if present
  if(!is.null(MLEobj$par)) {
    tmp = MARSScheckpar(MLEobj$par, n, m)
    if (!isTRUE(tmp)) msg = c(msg, tmp)
  }
  
  ## Check controls
  if(!is.null(MLEobj$control)){
    if(!is.list(MLEobj$control)) stop("Stopped in is.marssMLE() because control must be passed in as a list.\n", call.=FALSE)
  control = MLEobj$control
  en = names(alldefaults[[MLEobj$method]]$control)[!(names(alldefaults[[MLEobj$method]]$control) %in% c("boundsInits"))]

  for (el in en) {
    null.flag <- ( is.null(control[[el]]) && !(el %in% c("abstol") ) )  #abstol can be NULL
    if(null.flag) msg = c(msg, paste(el,"is missing from the control list\n"))

    if( !is.null(control[[el]]) ) {
      if( el %in% en[!(en %in% c("safe", "MCInit"))] ) {
        null.flag <- (!is.numeric(control[[el]]))
        if(null.flag) msg = c(msg, paste("control list element", el,"is non-numeric\n"))
      }
      if( (el %in% en[!(en %in% c("safe", "MCInit", "trace"))])  && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] <= 0)
        if(null.flag) msg = c(msg, paste("control list element", el,"less than or equal to zero\n"))
      }
      if (el %in% c("trace") && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] < 0)
        if(null.flag) msg = c(msg, paste("control list element", el,"less than zero\n"))
      }
      if (el %in% c("numInits", "numInitSteps", "trace", "minit", "maxit", "min.iter.conv.test", "conv.test.deltaT") && is.numeric(control[[el]])) {
        null.flag <- ( !is.wholenumber(control[[el]]) )
        if(null.flag) msg = c(msg, paste("control list element", el,"is not a whole number\n"))
      }
      if (el %in% c("minit") && !is.null(control$maxit) ) {
        if(!is.null(control$minit) && !is.null(control$maxit) && is.numeric(control$minit) && is.numeric(control$maxit) && is.wholenumber(control$minit) &&  is.wholenumber(control$maxit)) null.flag <- (control$minit > control$maxit)  
        if(null.flag) msg = c(msg, paste("control list element minit is greater than maxit\n"))
      }
      if (el %in% c("conv.test.deltaT") && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] < 2)
        if(null.flag) msg = c(msg, "control list element conv.test.deltaT must be greater than 1\n")
      }
      if (el %in% c("safe", "MCInit")) {
        null.flag <- !(control[[el]] %in% c(TRUE, FALSE) )	  
        if(null.flag) msg = c(msg, paste("control list element", el,"is not TRUE or FALSE\n"))
      }
    } # el is not null     
  }      # for el in en
  } #not null control
  
  ## Check control$boundsInits
    if(is.null(control$boundsInits)) msg = c(msg, "control$boundsInits is missing from the control list\n")
    if(!is.null(control$boundsInits)){
    en = c("B", "U", "logQ", "logR", "A", "Z")

    for (el in en) {
      target = control$boundsInits
      null.flag <- ( is.null(target[[el]]) )
      if(null.flag) msg = c(msg, paste("control$boundsInits list element", el,"is missing\n"))

      null.flag <- ( !is.null(target[[el]]) && length(target[[el]]) != 2)
      if(null.flag) msg = c(msg, paste("control$boundsInits list element", el,"is not a 2 element vector\n"))
 
      if (!null.flag) {
        null.flag <- (!is.numeric(target[[el]]))	  
        if(null.flag){ msg = c(msg, paste("control$boundsInits list element", el,"is not numeric\n"))
        }else {
           null.flag <- (target[[el]][1] >= target[[el]][2])
           if(null.flag) msg = c(msg, paste("The first element of control$boundsInits$", el," is not smaller than the second\n",sep=""))
           if (el == "B") {
              null.flag <- ( any(target[[el]] < 0) )	  
              if(null.flag) msg = c(msg, "One of the elements if control$boundsInits$B is < 0\n") }	  
        }
      }      
    }  
  } #if(!is.null(control$boundsInits))
    
  } # end if (isTRUE(msg))

if(length(msg) == 0){ return(TRUE)
}else {
  msg=c("\nErrors were caught in is.marssMLE()\n", msg)
  return(msg)
}
}
