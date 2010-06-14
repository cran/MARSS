#####
#
# is.marssMLE()
# Check that the marssMLE object has all the parts it needs 
# and that these have the proper size and form
#
#####

is.marssMLE <- function(MLEobj) 
{
  if(class(MLEobj) != "marssMLE") stop("is.marssMLE: Object class is not marssMLE.")
  msg = c()

  ## Check for required components
  el = c("model", "start", "control", "method")
  if( !all(el %in% names(MLEobj)) ) 
    msg = c(msg, paste("Element", el[!(el %in% names(MLEobj))], "is missing from object."))

  ## If required components present, check that the model is valid
  else msg = is.marssm(MLEobj$model)

  ## is.wholenumber() borrowed from is.integer example
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

  ## If model is OK, check inits and params
  if (isTRUE(msg)) {

  msg = list()

  ## Set m and n

  n = dim(MLEobj$model$fixed$Z)[1]
  m = dim(MLEobj$model$fixed$Z)[2]

  ## bug 584: check for T=1
  if( !is.null(dim(data)) && dim(data)[2] == 1 )
    msg = c(msg, "Data has only one time point.")

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
      msg = c(msg, paste("Missing or non-numeric initial value", en[init.null]))
    }
    if(any(dim.init)) {
      msg = c(msg, paste("Dimension problem: initial value", en[dim.init], "Dims should be", correct.dim1[dim.init], "x", correct.dim2[dim.init]))
    }    
  }

  ## Check params consistency if present

  if(!is.null(MLEobj$par)) {
    tmp = MARSScheckpar(MLEobj$par, n, m)
    if (!isTRUE(tmp)) msg = c(msg, tmp)
  }
  
  ## Check controls
  if(!is.null(MLEobj$control)){
    if(!is.list(MLEobj$control)) stop("is.marssMLE: control must be passed in as a list.")
  control = MLEobj$control
  control.null = NULL
  en = names(alldefaults[[MLEobj$method]]$control)[!(names(alldefaults[[MLEobj$method]]$control) %in% c("numInits", "numInitSteps", "boundsInits"))]

  for (el in en) {
    null.flag <- ( is.null(control[[el]]) )
 
    if (!null.flag) {
      if (el %in% c("abstol")) {
        null.flag <- (!is.numeric(control[[el]]) || control[[el]] <= 0)
      }
      if (el %in% c("iter.V0", "numInits", "numInitSteps", "trace")) {
        null.flag <- ( !is.wholenumber(control[[el]]) )
	null.flag <- ifelse(el == "trace", control[[el]] < 0, control[[el]] <= 0)
      }
      if (el %in% c("minit", "maxit")) {
        if(!is.wholenumber(control[[el]]) || control[[el]] <= 0){
          null.flag = TRUE
        }else if(!is.null(control$minit) && !is.null(control$maxit) && is.wholenumber(control$minit) &&  is.wholenumber(control$maxit)) {
	  null.flag <- (control$minit > control$maxit) }  
      }
      if (el %in% c("safe", "MCInit")) {
	null.flag <- !(control[[el]] %in% c(TRUE, FALSE) )	  
      }
    }      
    control.null <- c(control.null, null.flag)  
  }
  
  problem <- any(control.null)

  if (problem) {
    msg = c(msg, paste("Invalid control", en[control.null], "\n"))
  }
  }  

  ## Check control$boundsInits
  if(!is.null(control$boundsInits)){
    boundsInits.null = NULL
    en = c("B", "U", "logQ", "logR", "A", "Z")

    for (el in en) {
      target = control$boundsInits
      null.flag <- ( is.null(target[[el]]) || length(target[[el]]) != 2)
 
      if (!null.flag) {
        null.flag <- (!is.numeric(target[[el]]) || target[[el]][1] >= target[[el]][2])	  
        if (el == "B") {
          null.flag <- (null.flag || any(target[[el]] < 0))	  
        }
      }      
      boundsInits.null <- c(boundsInits.null, null.flag)  
    }  
    problem <- any(boundsInits.null)

    if (problem) {
      msg = c(msg, paste("Invalid control$boundsInits", en[boundsInits.null], "\n"))
    } 
  } #if(!is.null(control$boundsInits))
    
  } # end if (isTRUE(msg))

if(length(msg) == 0) return(TRUE)
else return(msg)
}
