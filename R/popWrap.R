#################################################################################
# popWrap() function
# Creates popWrap wrapper object (class popWrap). 
#################################################################################

## 'free' specifies which elements of a param are free and gives names for the free elements; elements that are equal
## will have the same name; elements in 'free' can be specified as numbers but they will be interpreted as names;
## non-free (i.e. fixed) values should be denoted with NA not 0 since the code will interpret 0 as the name "0" and assume
## that the user wants the params named "0" to be free.
## 'fixed' specifies the values of param elements that are fixed.  Values that are not fixed should be denoted NA.  
## All values in "fixed" must be numbers

# options for constraint
# "ignore" user is going to specify constraints with fixed and free
# Q & R  "unconstrained", "diagonal and unequal", "diagonal and all equal", "equalvarcov"
# U  "unconstrained", "all equal"
# B  "identity", "diagonal and unequal", "diagonal and all equal", "unconstrained"
# Z  "identity", a vector specifying which y goes to which x e.g. c(1, 1, 2, 1)
# A  "scaling"
# x0 "unconstrained" "equal" ("prior" not yet implemented)
# V0 No options, because it is determined by x0. if x0 is "fixed" then V0 must be 0; if x0 is prior then V0 must be prior)

# If constraint is passed in for R, B, U, Q or x0, then the corrsponding elements in fixed and free should not be passed in
# because 'constraint' is just a helper to set 'fixed' and 'free' for certain common cases.
# However, some params can be set with 'constraint' and others with 'fixed' and 'free'

popWrap <- function(y, allowed,
    inits=NULL,
    constraint=NULL,
    fixed=NULL, free=NULL, 
    miss.value=NULL,  
    control=NULL,
    method=NULL, 
    silent=FALSE)
{
  
### Specify some holders
wrapper.el = c("data", "m", "inits", "constraint", "fixed", "free", "miss.value", "control", "method")

req.args = c("inits", "constraint", "miss.value", "control", "method")
# model.elem   is specified in MARSSsettings

# Can't do anything if constraint is specified wrong
if(!is.null(constraint) && !is.list(constraint))
    stop("popWrap: constraint needs to be passed in as a list or left off (to use defaults)")

### If some elements are specified in both constraint and fixed or free, then use fixed/free and warn (b501)
both.cons.fixed.free = FALSE
for(i in unique(c(names(fixed), names(free))) )   {
      if(is.null(constraint)) constraint=list()
      if(!is.null(constraint[[i]])) both.cons.fixed.free = TRUE
      constraint[[i]] = "use fixed/free"
      } 
if(!silent && both.cons.fixed.free) warning("popWrap: both constraint and fixed (and/or free) specified for one of the model elements. fixed/free will be used and constraint arg ignored.")

defaults = alldefaults[[method]] #from MARSSsettings
## Now set defaults if needed, first deal with case where arg not passed in all
for(el in req.args)
  if(is.null(get(el))) assign(el, defaults[[el]])
	
### If some elements are missing from args use the defaults    
for(el in req.args) {
  tmp = get(el)
  if(!is.list(tmp) && !(el %in% c("miss.value", "method"))) stop(paste("popWrap: arg ",el," must be passed in as a list (or left off to use defaults)"))
  passed.in = (names(defaults[[el]]) %in% names(tmp))
  for(i in names(defaults[[el]])[!passed.in] )   {
    tmp[[i]] = defaults[[el]][[i]]
  }
  assign(el, tmp)
}

## KW Set m
n = ifelse(is.null(dim(y)), 1, dim(y)[1])
m = NA
if (identical(constraint$Z, "use fixed/free")) m = dim(fixed$Z)[2]
if (identical(constraint$Z, "identity")) m = n
if (is.factor(constraint$Z)) m = length(levels(constraint$Z)) 
if (is.matrix(constraint$Z)) m = dim(constraint$Z)[2] 

# if some of the elements of MCinit bounds were passed in but not others, use defaults
    default=defaults$control
    passed.in = names(default$boundsInits) %in% names(control$boundsInits)
    for(i in names(default$boundsInits)[!passed.in] )   {
      control$boundsInits[[i]] = default$boundsInits[[i]]
    }   

## Create popWrap object
  wrapperObj <- list(data=y, m=m, inits=inits, constraint=constraint, fixed=fixed, free=free, miss.value=miss.value, control=control, method=method)
  class(wrapperObj) <- "popWrap"

  ## checkPopWrap just checks that everything needed is present and
  ## the user didn't pass in any disallowed constraints that would prevent as.marss() from working
  ## mis-specified fixed/free will be caught by is.marssm()
  ## as.marssm() should construct the fixed/free from the constraint strings (unless "use fixed/free", then don't do anything)
  ## the exception is if either fixed or free is have no NAs, in that case, as.marssm() can construct the corresponding fixed/free
  ## e.g. fixed[[el]] = array(NA,dim=dim(fixed[[el]]))  Note dim not checked here, that will be caught by is.marssm()
  checkPopWrap(wrapperObj, wrapper.el, allowed, silent=silent)
  
  ## wrapperObj should now be ready for as.marssm()

  wrapperObj 
}
