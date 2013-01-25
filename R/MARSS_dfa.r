###################################################################################
# Helper function to create a DLM model sensu Zuur
# x(t)=x(t-1) + w(t), W~MVN(0,1)
# y(t)=Z x(t) + A(t) + D(t) d(t) + v(t), V~MVN(0,R)
# x(t0) = x0 + l, L ~ MVN(0,5)

# The conversion functions have 2 parts
# Part 1 Set up the DLM model in MARSS.marxss form
# Part 2 Call MARSS.marxss to finish the set-up and checking
###################################################################################
MARSS.dfa=function(MARSS.inputs){
# MARSS(data, model=list(), covariates=NULL, z.score=TRUE, demean=TRUE, control=list())
# model.defaults =list(A="zero", R="diagonal and equal", D="zero", x0="zero", V0=diag(5,1), tinitx=0, diffuse=FALSE, m=1)


#Part 1 Set up defaults and check that what the user passed in is allowed
# 1 Check for form dependent user inputs for method and reset defaults for inits, MCbounds, and control if desired
# 2 Specify the text shortcuts and whether factors or matrices can be passed in
#   The names in the allowed list do not need to be A, B, Q .... as used in the marssm object
#   Other names can be used if you want the user to use those names; then in the MARSS.form function
#   you convert the user passed in names into the marssm object with the A, B, Q, R, ... names
#   checkModelList() will check what the user passes in against these allowed values, so
#   so you need to make sure each name in model.defaults has a model.allowed value here 

#Set up some defaults
if(is.null(MARSS.inputs[["z.score"]])) MARSS.inputs$z.score=TRUE
if(is.null(MARSS.inputs[["demean"]])) MARSS.inputs$demean=TRUE

#Start error checking
problem=FALSE
msg=c()
#check that data and covariates elements are matrix or vector, no dataframes, and is numeric
for(el in c("data","covariates"[!is.null(MARSS.inputs[["covariates"]])])){
  if( !(is.matrix(MARSS.inputs[[el]]) || is.vector(MARSS.inputs[[el]])) ) {
    problem=TRUE
    msg = c(msg, paste(el," must be a matrix or vector (not a data frame or 3D array).\n"))
  }else{ 
    if(is.vector(MARSS.inputs[[el]])) MARSS.inputs[[el]]=matrix(MARSS.inputs[[el]],1)
  }
}
if(!is.null(MARSS.inputs$covariates)){
  if(dim(MARSS.inputs$data)[2]!= dim(MARSS.inputs$covariates)[2]){
    problem=TRUE
    msg = c(msg, "data and covariates must have the same number of time steps.\n")
  }
}
if(!is.na(MARSS.inputs$miss.value)){
  MARSS.inputs$data[MARSS.inputs$data==MARSS.inputs$miss.value]=NA
  if(!is.null(MARSS.inputs[["covariates"]])) MARSS.inputs$covariates[MARSS.inputs$covariates==MARSS.inputs$miss.value]=NA
  MARSS.inputs$miss.value=as.numeric(NA)  
}
if(!is.null(MARSS.inputs[["covariates"]])){
  if(any(is.na(MARSS.inputs$covariates))){
    problem=TRUE
    msg = c(msg, "covariates cannot have any missing values in a standard DFA.\n See User Guide section on DFA for alternate approaches when covariates have missing values.\n")
  }
}
if(!is.null(MARSS.inputs[["model"]][["m"]])){
  if(length(MARSS.inputs$model$m)!=1){
    problem=TRUE
    msg = c(msg, "model$m must be an integer between 1 and n.\n")
  }else{
  if( !is.numeric(MARSS.inputs$model$m) | !is.wholenumber(MARSS.inputs$model$m) ){
    problem=TRUE
    msg = c(msg, "model$m must be an integer between 1 and n.\n")
  }else{
  if(MARSS.inputs$model$m > dim(MARSS.inputs$data)[1]){
    problem=TRUE
    msg = c(msg, "model$m must be an integer between 1 and n.\n")
  }
  }
} }

if(problem){
    cat("\n","Errors were caught in MARSS.dfa \n", msg, sep="")
    stop("Stopped in MARSS.dfa() due to specification problem(s).\n", call.=FALSE)
}

n=dim(MARSS.inputs$data)[1]
model.allowed = list(
#if it is a length 1 vector then the value must be one of these.  All elements in your model list must be here
        A=c("unequal","zero"),
        R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        D=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        x0=c("unconstrained", "unequal", "zero"),
        V0=c("identity", "zero"),
        tinitx=c(0,1),
        diffuse=c(TRUE,FALSE),
        m=1:n,
#This line says what is allowed to be a matrix
        matrices = c("Z","A","R","D","x0","V0")
      )

if(is.null(MARSS.inputs[["model"]][["m"]])) m=1 else m=MARSS.inputs$model$m
if(!is.null(MARSS.inputs[["covariates"]])) D="unconstrained" else D="zero"

#defaults for any missing model list elements
model.defaults =list(A="zero", R="diagonal and equal", D=D, x0="zero", V0=diag(5,m), tinitx=0, diffuse=FALSE, m=1)

#This checks that what user passed in model list can be interpreted and converted to a marssm object
#if no errors, it updates the model list by filling in missing elements with the defaults
MARSS.inputs$model=checkModelList( MARSS.inputs$model, model.defaults, model.allowed )
model=MARSS.inputs$model

if(!(MARSS.inputs$z.score %in% c(TRUE,FALSE)))
   stop("Stopped in MARSS.dlm: demean must be TRUE/FALSE.\n", call.=FALSE)
if(!(MARSS.inputs$demean %in% c(TRUE,FALSE)))
   stop("Stopped in MARSS.dlm: demean must be TRUE/FALSE.\n", call.=FALSE)

#Set up Z
Z = matrix(list(), nrow=n, ncol=m)
# insert row (i) & col (j) indices
for(i in seq(n)) {Z[i,] = paste(i, seq(m), sep="")}
# set correct i,j values in Z to numeric 0
if(m>1) for(i in 1:(m-1)){Z[i,(i+1):m] = 0}

#Set up U
U=matrix(0,m,1)

#Set up D and d
if(is.null(MARSS.inputs[["covariates"]])) d=matrix(0,1,1) else d=MARSS.inputs$covariates

#Set up Q  & B
Q=diag(1,m); B=diag(1,m)

# set up list of model components for a marxss model
dfa.model = list(Z=Z, A=model$A, D=model$D, d=d, R=model$R, B=B, U=U, Q=Q, x0=model$x0, V0=model$V0)

dat=MARSS.inputs$data
if(MARSS.inputs$demean){
  y.bar = apply(dat, 1, mean, na.rm=TRUE)
  dat = (dat - y.bar)
}
if(MARSS.inputs$z.score){
  Sigma = sqrt(apply(dat, 1, var, na.rm=TRUE))
  dat = dat * (1/Sigma)
}

MARSS.inputs = list(data=dat, inits=MARSS.inputs$inits, MCbounds=MARSS.inputs$MCbounds, miss.value=as.numeric(NA), control=MARSS.inputs$control, method=MARSS.inputs$method, form="dlm", silent=MARSS.inputs$silent, fit=MARSS.inputs$fit)

tmp = MARSS.marxss(list(data=dat,model=dfa.model))
modelObj=tmp$marssm
modelObj$miss.value=as.numeric(NA)
MARSS.inputs$marssm = modelObj 
MARSS.inputs$form.info = tmp$form.info 

  ## Return MARSS inputs as list
MARSS.inputs
}

print_dfa = function(x){ return(print_marxss(x)) }
coef_dfa = function(x){ return(coef_marxss(x)) }
MARSSinits_dfa = function(MLEobj, inits){ return(MARSSinits_marxss(MLEobj, inits)) }
 