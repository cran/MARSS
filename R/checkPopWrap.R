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
  if( !all(el %in% names(wrapperObj)) ) 
    stop(paste("checkPopWrap: Element", el[!(el %in% names(wrapperObj))], "is missing from object."), call.=FALSE)

#Check constraints (b497)
  constraint = wrapperObj$constraint
  if (is.null(constraint)) { #previous line should have caught this
    stop("checkPopWrap: Constraint should not be NULL here. Look in popWrap for problem.", call.=FALSE)
  }
  else {
    if (!is.list(constraint)) stop("checkPopWrap: Constraint must be passed in as a list.", call.=FALSE)
    if( !all(names(constraint) %in% model.elem.w.V0 ) )
      stop("checkPopWrap: Incorrect element name(s) in constraint list (something misspelled?).", call.=FALSE)
    if( !all(model.elem %in% names(constraint)) ) 
      stop("checkPopWrap: constraint name missing.  This should have been caught in popWrap.", call.=FALSE)
  }
  
#check that constraint list doesn't have any duplicate names
  if(any(duplicated(names(constraint)))) stop("checkPopWrap: some of the constraint names are duplicated.", call.=FALSE)

#Series of checks on the model specification
problem = FALSE
msg=NULL
#check constraints only have allowed cases
  for (el in model.elem) {
    bad.str = FALSE
    #if length=1, then it must be a character string and that string must be in allowed. vectors length>1 are not allowed
    if(!is.factor(constraint[[el]]) && !is.matrix(constraint[[el]])) {
        if(length(constraint[[el]])!=1) bad.str=TRUE
        else if(!(constraint[[el]] %in% allowed[[el]]) ) bad.str=TRUE
    }
    if(bad.str) {
      problem=TRUE
      if(identical(constraint[[el]],"use fixed/free"))
      msg = c(msg, paste(" Use of fixed/free matrices to specify constraints is not allowed for ", el, "\n", sep=""))
      else msg = c(msg, paste(" The constraint value for ", el, " is not allowed.\n If you are trying to set a fixed parameter value, wrap value in matrix().\n If you are trying to set shared parameter elements, wrap in factor().\n If you are trying to set a text string, see ?MARSS for the allowed strings.\n", sep=""))
      }
    #if factor, must be allowed to be factor
    if(is.factor(wrapperObj$constraint[[el]]) && !(el %in% allowed$factors)) {
      problem=TRUE
      msg = c(msg, paste(" Constraint$",el," is not allowed to be a factor.\n", sep=""))
      }
    #if matrix must be allowed to be matrix
    if(is.matrix(wrapperObj$constraint[[el]]) && !(el %in% allowed$matrices)){
      problem=TRUE
      msg = c(msg, paste(" Constraint$",el," is not allowed to be a matrix.\n", sep=""))
      }     
  } # end for (el in model.elem)

  ## if fixed, free elements passed in, make sure they are matrix (otherwise as.marssm() will break) 
  for(mat in c("fixed","free")) {
    tmp=wrapperObj[[mat]]
    if(!is.null(tmp)) {
      if(!is.list(tmp)) {
        problem=TRUE
        msg = c(msg, paste(" ",mat," must be passed in as a list (or left off to ignore).\n"))
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
  if(identical(constraint$A,"scaling")){
   if(is.matrix(constraint$Z) && !is.design(constraint$Z)) {
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a design matrix:(0,1) and rowsums=1.\n")
        }
   if(identical(constraint$Z,"use fixed/free") || !is.null(wrapperObj$fixed$Z))
        if(!is.design(wrapperObj$fixed$Z)){
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a design matrix:(0,1) and rowsums=1.\n")
        }
   }
          
#check that data is numeric data
# bug 582
  if( !(is.matrix(wrapperObj$data)|| is.vector(wrapperObj$data)) ) {
    problem=TRUE
    msg = c(msg, " Data must be a matrix or vector (not a data frame).\n")
  }
  if(!is.numeric(wrapperObj$data)) {problem=TRUE; msg = c(msg, " Data must be numeric.\n")}

#check m; as.marssm() needs this
  if(!is.numeric(wrapperObj$m)) {
    problem=TRUE
    if(is.null(wrapperObj$m) || is.na(wrapperObj$m)) msg=c(msg, " m is not getting set.  Probably due to Z specification problem.\n")
    else msg=c(msg, " m (# state processes) must be numeric. m is set from the dim of Z.")
    }

#If there are errors, then don't proceed with the rest of the checks
  if(problem)  {
          cat("\n","Errors were caught in checkPopWrap \n", msg,"\n", sep="")
          stop("Stopped in checkPopWrap() due to specification problem(s).\n", call.=FALSE)
        }

#If a constraint is passed in as a factor, make sure it is the correct length otherwise as.marssm will break
problem = FALSE
msg=NULL
  n = ifelse(is.null(dim(wrapperObj$data)), 1, dim(wrapperObj$data)[1])
  m = wrapperObj$m

  ## Check constraints that are passed in as factors
  correct.factor.len = list(A=n, B=m, Q=m, R=n, U=m, x0=m, Z=n)  
  for (el in model.elem) {
    #if a factor then it needs to have the correct length otherwise as.marssm() won't work and no NAs
    if( is.factor(constraint[[el]]) ) {
      if(length(constraint[[el]]) != correct.factor.len[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" The constraint factor for ", el, " is not the right length. See help file.\n", sep=""))
          }
      if(NA.in.fac <- NA %in% constraint[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" NAs are not allowed in constraint factor for ", el, ". See help file.\n", sep=""))
          }
      }    
  } # end for (el in model.elem)
#if el == Z then factor needs to have m levels
    if( is.factor(constraint$Z) ) {
      if(length(levels(constraint$Z)) != m) {
            problem=TRUE
            msg=c(msg," When Z is a factor, the number of levels must equal the number of state processes (m).\n")
            }
    }  

#If there are errors in the factors
  if(problem)  {
          cat("\n","Errors were caught in checkPopWrap \n", msg,"\n", sep="")
          stop("Stopped in checkPopWrap() due to problem(s) with factors in constraint list.\n", call.=FALSE)
        }

  ## Check initial values consistency 
  inits = wrapperObj$inits
 
  # check that inits list doesn't have any duplicate names
  if(any(duplicated(names(inits)))) stop("Stopped in checkPopWrap: some of the inits names are duplicated.\n", call.=FALSE)

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
    cat("\n","Errors were caught in checkPopWrap \n", msg,"\n", sep="") 
    stop("Stopped in checkPopWrap() due to problem(s) with inits.\n", call.=FALSE)
  }
  
  return(TRUE)
}



