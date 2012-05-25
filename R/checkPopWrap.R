############################################################################################################################
#   checkPopWrap()
#   This checks user inputs to popWrap().
#   The main purpose is to make sure that as.marssm() will work not to make sure model is valid
#   Little error checking is done on fixed or free (e.g. no dim checks); that will be caught by is.marssm()
#   No error checking is done on controls and inits besides checking that it is present (NULL is ok); is.marssMLE() will error-check
##########################################################################################################################

checkPopWrap = function(wrapperObj, wrapper.el, allowed, silent=FALSE)
{
#####################################################################################
### Section 1, basic checks that nothing improper passed in
### and error in any will cause a stop because subsequent lines depend on previous lines
#####################################################################################

#Check that wrapper passed in and all wrapper elements are present
  if(!(class(wrapperObj)=="popWrap")) stop("checkPopWrap: Object class is not popWrap.", call.=FALSE)
  el = wrapper.el
  if( !all(el %in% names(wrapperObj)) ){ 
     msg=paste(" Element", el[!(el %in% names(wrapperObj))], "is missing from object.\n")
     cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
     stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
    }
#Check model structures (b497)
  model = wrapperObj$model
  if (is.null(model)){ #previous line should have caught this
     cat("\n","Errors were caught in checkPopWrap \n", " model should not be NULL here. Look in popWrap for problem.\n", sep="")
     stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
  }
  else {
    if (!is.list(model)){
      msg=" model must be passed in as a list.\n"
      cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
      stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
    }
    if( !all(names(model) %in% model.elem ) ){
      msg=" Incorrect element name(s) in model list (something misspelled?).\n"
      cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
      stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
      }
    if( !all(model.elem %in% names(model)) ){ 
      msg=" A model name is missing.  This should have been caught in popWrap.\n"
      cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
      stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
      }
  }
  
#check that model list doesn't have any duplicate names
  if(any(duplicated(names(model)))){
      msg=" Some of the model names are duplicated.\n"
      cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
      stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
   }

#Series of checks on the model specification
problem = FALSE
msg=NULL
#check model structures only have allowed cases
  for (el in model.elem) {
    bad.str = FALSE
    #if length=1, then it must be a character string and that string must be in allowed. vectors length>1 are not allowed
    if(!is.factor(model[[el]]) && !is.matrix(model[[el]])) {
        if(length(model[[el]])!=1) bad.str=TRUE
        else if(!(model[[el]] %in% allowed[[el]]) ) bad.str=TRUE
    }
    if(bad.str) {
      problem=TRUE
      msg = c(msg, paste(" The model value for ", el, " is not allowed.\n See ?MARSS for the allowed values.\n", sep=""))
      }
    #if factor, must be allowed to be factor
    if(is.factor(wrapperObj$model[[el]]) && !(el %in% allowed$factors)) {
      problem=TRUE
      msg = c(msg, paste(" model$",el," is not allowed to be a factor.\n", sep=""))
      }
    #if matrix must be allowed to be matrix
    if(is.matrix(wrapperObj$model[[el]]) && !(el %in% allowed$matrices)){
      problem=TRUE
      msg = c(msg, paste(" model$",el," is not allowed to be a matrix.\n", sep=""))
      }     
    #if matrix then no NAs if character or list; this would be caught in is.marssm but would be hard for user to understand problem
    if(is.matrix(wrapperObj$model[[el]]) && (el %in% allowed$matrices)){
      if( is.character(wrapperObj$model[[el]]) && any(is.na(wrapperObj$model[[el]])) ){
        problem=TRUE
        msg = c(msg, paste(" model$",el," is a character matrix. No NAs allowed in this case.\n", sep=""))
        }
      if( is.list(wrapperObj$model[[el]]) && any(is.na(wrapperObj$model[[el]])) ){
        problem=TRUE
        msg = c(msg, paste(" model$",el," is a list matrix. No NAs allowed in this case.\n", sep=""))
        }
      if( is.list(wrapperObj$model[[el]]) && any(sapply(wrapperObj$model[[el]],length)!=1) ){
        problem=TRUE
        msg = c(msg, paste(" model$",el," is a list matrix. No NAs allowed in this case.\n", sep=""))
        }

      } # is matrix  
    } # end for (el in model.elem)

  ## if fixed, free elements passed in, make sure they are matrix (otherwise as.marssm() will break) 
  for(mat in c("fixed","free")) {
    tmp=wrapperObj[[mat]]
    if(!is.null(tmp)) {
      if(!is.list(tmp)) {
        problem=TRUE
        msg = c(msg, paste(" ",mat," must be passed in as a list (or left off to ignore).\n",sep=""))
        }
      passed.in = (model.elem %in% names(tmp))
      for(i in model.elem[passed.in])
        if(!is.matrix(tmp[[i]])) {
        problem = TRUE
        msg = c(msg, paste(" ",mat,"$",i," must be passed in as a matrix.\n"))
        }
    }
  }  #check of fixed and free
  
  #Check that if A is scaling, then Z spec will lead to a design matrix b589
  if(identical(model$A,"scaling")){
   if(is.matrix(model$Z) && !is.design(model$Z)) {
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a design matrix:(0,1) and rowsums=1.\n")
        }
   }
          
#check that data is numeric data
# bug 582
  if( !(is.matrix(wrapperObj$data) || is.vector(wrapperObj$data)) ) {
    problem=TRUE
    msg = c(msg, " Data must be a matrix or vector (not a data frame).\n")
  }
  if(!is.numeric(wrapperObj$data)) {problem=TRUE; msg = c(msg, " Data must be numeric.\n")}

#check that miss.value is numeric or NA; make sure R doesn't think miss.value=NA is logical
  if(is.na(wrapperObj$miss.value)) wrapperObj$miss.value = as.numeric(wrapperObj$miss.value)
  if( !(length(wrapperObj$miss.value)==1) || !(is.na(wrapperObj$miss.value) || is.numeric(wrapperObj$miss.value))){
    problem=TRUE
    msg = c(msg, " miss.value must be length 1 and numeric (or NA).\n")
  }

#check m; as.marssm() needs this
  if(!is.numeric(wrapperObj$m)) {
    problem=TRUE
    if(is.null(wrapperObj$m) || is.na(wrapperObj$m)) msg=c(msg, " m is not getting set.  Probably due to Z specification problem.\n")
    else msg=c(msg, " m (# state processes) must be numeric. m is set from the dim of Z.")
    }

#If there are errors, then don't proceed with the rest of the checks
  if(problem)  {
          cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
          stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
        }

#If a model is passed in as a factor, make sure it is the correct length otherwise as.marssm will break
problem = FALSE
msg=NULL
  n = ifelse(is.null(dim(wrapperObj$data)), 1, dim(wrapperObj$data)[1])
  m = wrapperObj$m

  ## Check model structures that are passed in as factors
  correct.factor.len = list(A=n, B=m, Q=m, R=n, U=m, x0=m, Z=n)  
  for (el in model.elem) {
    #if a factor then it needs to have the correct length otherwise as.marssm() won't work and no NAs
    if( is.factor(model[[el]]) ) {
      if(length(model[[el]]) != correct.factor.len[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" The model$", el, " is being passed as a factor and should be length ", correct.factor.len[[el]], " based on data dims and Z. It's not. See help file.\n", sep=""))
          }
      if(NA.in.fac <- NA %in% model[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" NAs are not allowed in model factor for ", el, ". See help file.\n", sep=""))
          }
      }    
  } # end for (el in model.elem)
#if el == Z then factor needs to have m levels
    if( is.factor(model$Z) ) {
      if(length(levels(model$Z)) != m) {
            problem=TRUE
            msg=c(msg," When Z is a factor, the number of levels must equal the number of state processes (m).\n")
            }
    }  

#If there are errors in the factors
  if(problem)  {
          cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
          stop("Stopped in checkPopWrap() due to problem(s) with factors in model list.\n", call.=FALSE)
        }

  ## Check initial values consistency 
  inits = wrapperObj$inits
 
  # check that inits list doesn't have any duplicate names
  if(any(duplicated(names(inits)))){
     msg = " Some of the inits names are duplicated.\n"
     cat("\n","Errors were caught in checkPopWrap \n", msg, sep="")
     stop("Stopped in checkPopWrap() due to problem(s) with inits.\n", call.=FALSE)
  }

  problem = FALSE
  msg=NULL
  init.null = dim.init = NULL
  # Z is dealt with later
  elem = model.elem[model.elem != "Z"]
  
  for (el in elem) {
    init.null.flag <- ( is.null(inits[[el]]) || !is.numeric(inits[[el]]) )
    dim.init.flag = FALSE
 
    if (!init.null.flag) {  
      if(el %in% c("B", "Q", "V0")) {
        if(is.matrix(inits[[el]])) dim.init.flag <- !( dim(inits[[el]])[1]==m && dim(inits[[el]])[2]==m )
        else dim.init.flag <- !(length(inits[[el]])==1)
      }
      if(el %in% c("U", "x0")) {
	if(is.matrix(inits[[el]])) dim.init.flag <- !( dim(inits[[el]])[1]==m && dim(inits[[el]])[2]==1 )
        else dim.init.flag <- !(length(inits[[el]])==1)
      }
      if(el == "A") {
	if(is.matrix(inits[[el]])) dim.init.flag <- !( dim(inits[[el]])[1]==n && dim(inits[[el]])[2]==1 )
        else dim.init.flag <- !(length(inits[[el]])==1)
      }
      if(el == "R") {
	if(is.matrix(inits[[el]])) dim.init.flag <- !( dim(inits[[el]])[1]==n && dim(inits[[el]])[2]==n )
        else dim.init.flag <- !(length(inits[[el]])==1)
      }
    }
 
    init.null <- c(init.null, init.null.flag)
    dim.init <- c(dim.init, dim.init.flag)  
  }
  
  problem <- any(c(init.null, dim.init))

  if (problem) {    
    if(any(init.null)) {
      msg = c(msg, paste("Missing or non-numeric initial value", elem[init.null], "\n"))
    }
    if(any(dim.init)) {
      msg = c(msg, paste("Dimension problem: initial value", elem[dim.init], "\n"))
    }
    cat("\n","Errors were caught in checkPopWrap \n", msg, sep="") 
    stop("Stopped in checkPopWrap() due to problem(s) with inits.\n", call.=FALSE)
  }
  
  return(TRUE)
}
