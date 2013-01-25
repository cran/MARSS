###################################################################################
# Coerce model list input in MARSS() call to MARSS model object (class marssm)
# using form = MARXSS
# x(t)=B(t)x(t-1) + U(t) + C(t)c(t) + w(t), W~MVN(0,Q)
# y(t)=Z(t)x(t) + A(t) + D(t)d(t) + v(t), V~MVN(0,R)
# x(t0) = x0 + l, L ~ MVN(0,V0)
# produces model object with fixed, free, miss.value, X.names, tinitx and data

# The conversion functions have 2 parts
# Part 1 Set up the defaults and allowed structures
# Part 2 Do the conversion of model list to a marssm object
###################################################################################
MARSS.marxss=function(MARSS.call){
#Part 1 Set up defaults and check that what the user passed in is allowed
# This part is not required.  All that is required is that a proper marssm object is returned

# 1 Check for form dependent user inputs for method and reset defaults for inits, MCbounds, and control if desired

# 2 Specify the text shortcuts and whether factors or matrices can be passed in
#   The names in the allowed list do not need to be A, B, Q .... as used in the marssm object
#   Other names can be used if you want the user to use those names; then in the MARSS.form function
#   you convert the user passed in names into the marssm object with the A, B, Q, R, ... names
#   checkModelList() will check what the user passes in against these allowed values, so
#   so you need to make sure each name in model.defaults has a model.allowed value here 
model.allowed = list(
        A=c("scaling", "unconstrained", "unequal", "equal", "zero"),
        D=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal","zero"),
        B=c("identity", "zero", "unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Q=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        U=c("unconstrained", "equal", "unequal", "zero"),
        C=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal","zero"),
        x0=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal"),
        V0=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Z=c("identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "onestate"),
        c=c("zero"),
        d=c("zero"),
        tinitx=c(0,1),
        diffuse=c(TRUE,FALSE),
        factors = c("Z"),
        matrices = c("A","B","Q","R","U","x0","Z","V0","D","C","d","c")
      )

#model.defaults is form dependent so you must specify it
model.defaults =list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero", D="zero", d=matrix(0,1,1), C="zero", c=matrix(0,1,1), tinitx=0, diffuse=FALSE)

if(!is.null(MARSS.call[["model"]][["c"]]))
  if(!identical(MARSS.call$model$c, "zero") & !all(MARSS.call$model$c==0)) model.defaults$C="unconstrained"
if(!is.null(MARSS.call[["model"]][["d"]]))
  if(!identical(MARSS.call$model$d, "zero") & !all(MARSS.call$model$d==0)) model.defaults$D="unconstrained"

#This checks that what user passed in model list can be interpreted and converted to a marssm object
#if no errors, it updates the model list by filling in missing elements with the defaults
MARSS.call$model=checkModelList( MARSS.call$model, model.defaults, model.allowed )

# Part 2 Convert the model list to a marssm object  
  ## set up fixed and free elements
  fixed = free = list()

  model = MARSS.call$model
  model.elem = c("Z","A","R","B","U","Q","x0","V0","D","d","C","c") 
  dat = MARSS.call$data 
  if(is.vector(dat)) dat=matrix(dat,1,length(dat))
  if(is.null(model[["X.names"]]) & identical(model$Z,"identity"))
    model$X.names=rownames(dat)
  if(is.null(model[["X.names"]]) & is.matrix(model$Z))
    if(is.identity(model$Z)) model$X.names=rownames(dat)
  n = dim(dat)[1]; TT = dim(dat)[2]
  ## Set m based on Z specification IF Z was specified; errors will be reported later if m conflicts with other parameters
  m = NA
  if (identical(model$Z, "unconstrained")) m = n
  if (identical(model$Z, "equalvarcov")) m = n
  if (identical(model$Z, "diagonal and equal")) m = n
  if (identical(model$Z, "diagonal and unequal")) m = n
  if (identical(model$Z, "onestate")) m = 1
  if (identical(model$Z, "identity")) m = n
  if (is.factor(model$Z)) m = length(levels(model$Z)) 
  if (is.array(model$Z)) m = dim(model$Z)[2] 

  ## Set c1 based on C specification if C specified with a particular shape
  ## error checking later will complain if C and c (or D and d) conflict
  if (is.array(model$C)) c1 = dim(model$C)[2] else c1 = 1
  if (is.array(model$D)) d1 = dim(model$D)[2] else d1 = 1

  for(el in c("c","d")){
    thedim=get(paste(el,"1",sep=""))
    if(identical(model[[el]], "zero")){ model[[el]]=matrix(0,thedim,1) }
    #if(!any(is.na(model[[el]])) & all(model[[el]]==0)) model[[toupper(el)]]="zero"
    if(is.vector(model[[el]])) model[[el]]=matrix(model[[el]],1,length(model[[el]]))
  }
  #Now set c1 and d1 based on c and d, which should now be a matrix of some sort.  This ensures that c1 and d1 are set
  c1=dim(model$c)[1]; d1=dim(model$d)[1]
  model.dims = list(data=c(n,TT),Z=c(n,m),U=c(m,1),A=c(n,1),B=c(m,m),Q=c(m,m),R=c(n,n),x0=c(m,1),V0=c(m,m),D=c(n,d1), d=c(d1,1), C=c(m,c1), c=c(c1,1))

  ## Error-checking section that is specific to marxss form
  # Note most error checking happens in checkMARSSInputs, checkModelList, and is.marssMLE
  #If a model is passed in as a factor, make sure it is the correct length otherwise construction of marssm object will break
  problem = FALSE
  msg=NULL
  ## Check model structures that are passed in as factors
  correct.factor.len = list(Z=n,A=n,R=n,B=m,U=m,Q=m,x0=m,V0=m)  
  for (el in model.elem) {
    #if a factor then it needs to have the correct length otherwise construction of marssm object will break and no NAs
    if( is.factor(model[[el]]) ) {
      if(length(model[[el]]) != correct.factor.len[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" The model$", el, " is being passed as a factor and should be length ", correct.factor.len[[el]], " based on data dims. It's not. See help file.\n", sep="")) }
      if(NA.in.fac <- NA %in% model[[el]]) {
          problem=TRUE
          msg = c(msg, paste(" NAs are not allowed in model factor for ", el, ". See help file.\n", sep="")) }
      }  #is factor  
  } # end for (el in model.elem)

#if el == Z then factor needs to have m levels
    if( is.factor(model$Z) ) {
      if(length(levels(model$Z)) != m) {
            problem=TRUE
            msg=c(msg," When Z is a factor, the number of levels must equal the number of state processes (m).\n")
      } }

#Check that if A is scaling, then Z spec must lead to a design matrix
  if(identical(model$A,"scaling")){
   if((!is.array(model$Z) & !is.factor(model$Z))) #if it is a string
     if(!(model$Z %in% c("onestate","identity"))){
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a design matrix:(0,1) and rowsums=1.\n")
     }
   if(is.matrix(model$Z) && !is.design(model$Z)) { #if it is a matrix
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\n")
        }
   if(is.array(model$Z))
     if(dim(model$Z)[3]!=1 && !is.design(model$Z)) { #if it is a matrix
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\n")
        }
   }  
   if(is.array(model$x0) & length(dim(model$x0))==3)
      if(dim(model$x0)[3]!=1){
        problem=TRUE
        msg = c(msg, " x0 cannot be time-varying thus if x0 in model arg is 3D, the 3rd dim must equal 1.\n")
      } 
   if(is.array(model$V0) & length(dim(model$V0))==3)
      if(dim(model$V0)[3]!=1){
        problem=TRUE
        msg = c(msg, " V0 cannot be time-varying thus if V0 in model arg is 3D, the 3rd dim must equal 1.\n")
      } 
 #if C is diagonal and equal or diagonal and unequal, then d1=m
   if(identical(model$C,"diagonal and equal") | identical(model$C,"diagonal and unequal"))
     if(c1!=m){
        problem=TRUE
        msg = c(msg, " If C is diagonal, it must be square and c must be m x 1.\n")
     }
 #if D is diagonal and equal or diagonal and unequal, then d1=n
   if(identical(model$D,"diagonal and equal") | identical(model$D,"diagonal and unequal"))
     if(d1!=n){
        problem=TRUE
        msg = c(msg, " If D is diagonal, it must be square and d must be n x 1.\n")
     }
 #if c and d can't have any NAs or Infs
   for(el in c("c","d"))
     if(any(is.na(model[[el]])) | !is.numeric(model[[el]]) | any(is.infinite(model[[el]]))){
        problem=TRUE
        msg = c(msg, paste(" ",el,"must be numeric and have no NAs, NaNs, or Infs.\n"))
     }
 #c and d must be a 2D matrix and 2nd dim must be 1 or TT
   for(el in c("c","d")){
     if( length(dim(model[[el]]))!=2 ){
        problem=TRUE
        msg = c(msg, paste(" ",el,"must be a 2D matrix.\n"))
     }else{
       if( !(dim(model[[el]])[2] == 1 | dim(model[[el]])[2] == TT) ){
         problem=TRUE
         msg = c(msg, paste(" ",el,"must be a 2D matrix with 2nd dim equal to 1 or T (length of data).\n"))
       }
     }
   }
          
#If there are problems
  if(problem)  {
          cat("\n","Errors were caught in MARSS.marxss \n", msg, sep="")
          stop("Stopped in MARSS.marxss() due to problem(s) with model specification.\n", call.=FALSE)
        }
  #end of error section       

  ##Translate the wrapper object into the marssm object 
  ## Translate the model structure names (shortcuts) into fixed and free
  ## fixed is a dim(1)*dim(2) X 1 vector of the fixed (intercepts) values
  ## free is a dim(1)*dim(2) X p vector of the free (betas) values for the p estimated elements
  
  model.elem = c("Z","A","R","B","U","Q","x0","V0","D","C")
  if(which(model.elem=="Z")>which(model.elem=="A")) model.elem=rev(model.elem) #Z must go first

  tmp=list()
  for(el in model.elem) {
   tmp[[el]]="not assigned"
   if(el=="Z" & is.factor(model$Z)) { 
      X.names = unique(model$Z)
      tmp[[el]] = matrix(0,model.dims$Z[1], model.dims$Z[2])  
      for(i in X.names) tmp[[el]][which(model$Z==i), which(as.vector(X.names)==i)] = 1
    }
    if(el=="Z" & identical(model$Z,"onestate")) {   #m=1
      tmp[[el]] = matrix(1, n, 0)
    }
   if( identical(model[[el]],"identity") ) { 
    tmp[[el]] = diag(1,model.dims[[el]][1])
    }
   if(identical(model[[el]],"diagonal and equal")) {
    tmp[[el]] = array(list(0),dim=model.dims[[el]])
    diag(tmp[[el]])="diag" #paste(el,"(diag)",sep="")
    if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
   }
   if(identical(model[[el]],"diagonal and unequal")) {
    tmp[[el]] = array(list(0),dim=model.dims[[el]])
    dim.mat = model.dims[[el]][1]
    diag(tmp[[el]])=paste("(",as.character(1:dim.mat),",",as.character(1:dim.mat),")",sep="") #paste(el,"(",as.character(1:dim.mat),",",as.character(1:dim.mat),")",sep="")
    if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
  }
  if(identical(model[[el]],"unconstrained") | identical(model[[el]],"unequal")){
    tmp[[el]]=array(NA,dim=model.dims[[el]]) 
    if(el %in% c("Q","R","V0")){  #variance-covariance matrices
      dim.mat = model.dims[[el]][1]
      for(i in 1:dim.mat){
        for(j in 1:dim.mat) tmp[[el]][i,j]=tmp[[el]][j,i]=paste("(",i,",",j,")",sep="") #paste(el,"(",i,",",j,")",sep="")
      }
    }else{ #not var-cov matrix
      row.name=1:model.dims[[el]][1]
      col.name=1:model.dims[[el]][2]
      if(el %in% c("C","D")){
        if(el=="C" & !is.null(model[["X.names"]])) row.name=model$X.names
        if(el=="D" & !is.null(rownames(dat))) row.name=rownames(dat)
        if(!is.null(rownames(model[[tolower(el)]]))) col.name=rownames(model[[tolower(el)]])
      }
      for(i in 1:model.dims[[el]][1]){
        for(j in 1:model.dims[[el]][2]){
         if(model.dims[[el]][2]>1) tmp[[el]][i,j]=paste("(",row.name[i],",",col.name[j],")",sep="") #paste(el,"(",row.name[i],",",col.name[j],")",sep="")
         else tmp[[el]][i,j]=paste(row.name[i],sep=",") #paste(el,row.name[i],sep=",")
         }
      }
    }
    if(length(tmp[[el]])==1) tmp[[el]][1,1]=el 
  } #unconstrained
  if(identical(model[[el]],"equalvarcov")) {
    tmp[[el]]=array("offdiag",dim=model.dims[[el]]) #array(paste(el,"(offdiag)",sep=""),dim=model.dims[[el]])
    diag(tmp[[el]])="diag" #paste(el,"(diag)",sep="")
    if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
  }
  if(identical(model[[el]],"equal")) { 
    tmp[[el]]=array("1",dim=model.dims[[el]]) #array(el,dim=model.dims[[el]])
  }
  if(identical(model[[el]],"zero")) { 
    tmp[[el]]=array(0,dim=model.dims[[el]])
  }
  if(is.array(model[[el]])) {
    tmp[[el]]=model[[el]]
    }
  if(el=="A" & identical(model[[el]], "scaling")) {  #check above ensures that Z is design and time-invariant
    ## Construct A from fixed Z matrix
    tmp[[el]] = matrix(list(),model.dims$A[1],model.dims$A[2])
    tmp[[el]][,1]=as.character(1:model.dims$A[1])  
    for(i in 1:m) {
      tmp[[el]][min(which(tmp$Z[,i]!=0)), 1] = 0
    }
  }     
  if(identical(tmp[[el]],"not assigned")) stop(paste("MARSS.marxss: tmp was not assigned for ",el,".\n",sep=""))
  free[[el]] = convert.model.mat(tmp[[el]])$free
  fixed[[el]] = convert.model.mat(tmp[[el]])$fixed    
}
marxss_object = list(fixed=fixed, free=free, data=dat, miss.value=MARSS.call$miss.value, X.names=model$X.names, tinitx=model$tinitx, diffuse=model$diffuse, model.dims=model.dims, form="marxss" )

#This is the f+Dp form for the MARXSS model used for user displays, printing and such
#Below this gets converted to the base MARSS form used in the algorithms
MARSS.call$form.info$marxss = list(free=free, fixed=fixed)
MARSS.call$form.info$model.dims = model.dims #because printing will need to know the dims of the D, C, d, and c matrices
MARSS.call$form.info$form = "marxss" #This is not what MARSS() was called with but rather what form.info is based on

#This step converts U+Cc into equivalent Uu and A+Dd into Aa
#So U --> [C U] and u --> [c \\ 1]; A --> [D A] and a --> [d \\ 1]
  for(el in c("C","D")){
    if(dim(tmp[[el]])[1]!=model.dims[[el]][1] | dim(tmp[[el]])[2]!=model.dims[[el]][2])
      stop("Stopped in MARSS.marxss(): ",el," should have ",model.dims[[el]][1], " x ", model.dims[[el]][2]," dims and it doesn't.\n", call.=FALSE)
    if(dim(model[[tolower(el)]])[1]!=model.dims[[tolower(el)]][1] )
      stop("Stopped in MARSS.marxss(): ",tolower(el)," should have ",model.dims[[tolower(el)]][1], " rows and it doesn't.\n", call.=FALSE)
    if(el=="C") el2="U" else el2="A"
    if(dim(tmp[[el2]])[1]!=model.dims[[el2]][1] | dim(tmp[[el2]])[2]!=model.dims[[el2]][2])
      stop("Stopped in MARSS.marxss(): ",el2," should have ",model.dims[[el2]][1], " x ", model.dims[[el2]][2]," dims and it doesn't.\n", call.=FALSE)
    if(!all(unlist(lapply(unname(tmp[[el]]),identical, 0)))){ #works if C/D is matrix(0,m,0) too
      model[[tolower(el2)]]=rbind(model[[tolower(el)]],1) #so put in u or a
      Tmax.fixed=max(dim(fixed[[el]])[3],dim(fixed[[el2]])[3])
      Tmax.free=max(dim(free[[el]])[3],dim(free[[el2]])[3])
      dim.fixed=c(dim(fixed[[el]])[1],dim(fixed[[el]])[2],dim(fixed[[el2]])[1],dim(fixed[[el2]])[2])
      dim.free=c(dim(free[[el]])[1],dim(free[[el]])[2],dim(free[[el2]])[1],dim(free[[el2]])[2])
      tmp.fixed=array(NA,dim=c(dim.fixed[1]+dim.fixed[3],1,Tmax.fixed))
      tmp.free=array(0,dim=c(dim.free[1]+dim.free[3],dim.free[2]+dim.free[4],Tmax.free))
      tmp.fixed[1:dim.fixed[1],,]=fixed[[el]]
      tmp.fixed[(dim.fixed[1]+1):dim(tmp.fixed)[1],,]=fixed[[el2]]
      if(dim.free[2]>0) tmp.free[1:dim.free[1],1:dim.free[2],]=free[[el]]
      if(dim.free[4]>0) tmp.free[(dim.free[1]+1):dim(tmp.free)[1],(dim.free[2]+1):dim(tmp.free)[2],]=free[[el2]]
      colnames(tmp.free)=c(colnames(free[[el]]),colnames(free[[el2]]))
      fixed[[el2]]=tmp.fixed
      free[[el2]]=tmp.free
      model.dims[[el2]][2]=model.dims[[el]][2]+model.dims[[el2]][2]
    }else{ model[[tolower(el2)]]=matrix(1,1,1) }
  }

#this part converts U(t)u(t) to U(t) and A(t)a(t) to A(t); the small case are inputs and the large case are estimated parameters  
  for(el in c("U","A")){
  if(!identical(unname(model[[tolower(el)]]), matrix(1,1,1))){
    if(dim(free[[el]])[1]!=model.dims[[el]][1]*model.dims[[el]][2])
       stop("Stopped in MARSS.marxss(): free$",el," should have ",model.dims[[el]][1]*model.dims[[el]][2]," rows and it doesn't.\n", call.=FALSE)
    if(dim(fixed[[el]])[1]!=model.dims[[el]][1]*model.dims[[el]][2])
       stop("Stopped in MARSS.marxss(): fixed$",el," should have ",model.dims[[el]][1]*model.dims[[el]][2]," rows and it doesn't.\n", call.=FALSE)
    free.orig=free[[el]]; fixed.orig=fixed[[el]]
    dim.free2=dim(free.orig)[2]; dim.free3=dim(free.orig)[3];  dim.fixed3=dim(fixed.orig)[3]; 
    dim.u.3=dim(model[[tolower(el)]])[2]
    Tmax=max(dim.free3, dim.fixed3, dim.u.3)
    free[[el]]=array(0,dim=c(model.dims[[el]][1],dim.free2,Tmax))
    colnames(free[[el]])=colnames(free.orig)
    fixed[[el]]=array(0,dim=c(model.dims[[el]][1],1,Tmax))
    for(t in 1:Tmax){
      f.t=sub3D(fixed.orig,t=min(t,dim.fixed3))
      d.t=sub3D(free.orig,t=min(t,dim.free3))
      ua.t=model[[tolower(el)]][,min(t,dim.u.3),drop=FALSE]      
      free[[el]][,,t]=(t(ua.t)%x%diag(1,model.dims[[el]][1]))%*%d.t
      fixed[[el]][,,t]=(t(ua.t)%x%diag(1,model.dims[[el]][1]))%*%f.t
    }
  }
  }
  marssm.elem = c("Z","A","R","B","U","Q","x0","V0")
  marssm.dims = list(data=c(n,TT),Z=c(n,m),U=c(m,1),A=c(n,1),B=c(m,m),Q=c(m,m),R=c(n,n),x0=c(m,1),V0=c(m,m))
free=free[marssm.elem]
  fixed=fixed[marssm.elem]
  
  #Save the X names coming in from model$Z otherwise this information is lost
  if( is.null(model[["X.names"]]) ){
    model$X.names=paste("X",1:m,sep="")
    if(is.array(model$Z) & !is.null(colnames(model$Z))) model$X.names=colnames(model$Z)
    if(is.factor(model$Z)) model$X.names=unique(model$Z)
  } 
  
  ## Create marssm object
  modelObj = list(fixed=fixed, free=free, data=dat, miss.value=MARSS.call$miss.value, X.names=model$X.names, tinitx=model$tinitx, diffuse=model$diffuse, model.dims=marssm.dims, form="marssm" )
  class(modelObj) = "marssm"
  MARSS.call$marssm=modelObj
  ## Return MARSS call list with marssm obj added
  MARSS.call
}

marssm_to_marxss=function(x){
  #This changes the model part of the MLE object for printing purposes.  It will act as if the user specified a 
  #MARXSS instead of MARSS
  x$model$fixed = x$form.info$marxss$fixed
  x$model$free = x$form.info$marxss$free
  for(val in c("par","start","par.se","par.bias","par.upCI","par.lowCI")){
  if(!is.null(x[[val]])){
  tmp.dim=dim(x$form.info$marxss$free$C)[2]
  if(tmp.dim==0){
   x[[val]][["C"]] = matrix(0,0,1)
  }else{
   #because marssm.U is [marxss.C marxss.U]
   x[[val]][["C"]] = x[[val]][["U"]][1:tmp.dim,, drop=FALSE]
   x[[val]][["U"]] = x[[val]][["U"]][-(1:tmp.dim),, drop=FALSE]
  }
  tmp.dim=dim(x$form.info$marxss$free$D)[2]
  if(tmp.dim==0){
   x[[val]][["D"]] = matrix(0,0,1)
  }else{
   #because marssm.A is [marxss.D marxss.A]
   x[[val]][["D"]] = x[[val]][["A"]][1:tmp.dim,, drop=FALSE]
   x[[val]][["A"]] = x[[val]][["A"]][-(1:tmp.dim),, drop=FALSE]
  }
  }
  } #for val in par, start
  x$model$d = x$call$model$d
  x$model$c = x$call$model$c
  x$model$model.dims = x$form.info$model.dims
  x$model$form = "marxss"
  return(x)
}

print_marxss = function(x){ return(marssm_to_marxss(x)) }

coef_marxss = function(x){ return(marssm_to_marxss(x)) }

MARSSinits_marxss = function(MLEobj, inits){
  #B, Z, R, Q, x0 and V0 stay the same
  #U and A change
  #this function will return a U and A element for inits
  if(is.null(inits)) inits=list()
  inits.defaults=list(
    U=alldefaults[[MLEobj$method]][["inits"]][["U"]], C=0,
    A=alldefaults[[MLEobj$method]][["inits"]][["A"]], D=0 )
  elems=c("U","A","C","D")
  for(elem in elems){
    tmp.dim=dim(MLEobj$form.info$marxss$free[[elem]])[2]
    if(!is.null(inits[[elem]])){
      if(!(length(inits[[elem]]) %in% c(tmp.dim,1))) 
         stop(paste("MARSSinits: ", elem," inits must be either a scalar (dim=NULL) or a matrix with 1 col and rows equal to the num of est values in ",elem,".",sep=""), call.=FALSE )
      if(tmp.dim!=0) inits[[elem]] = matrix(inits[[elem]],tmp.dim,1) else inits[[elem]]=matrix(0,0,1)
    }else{  inits[[elem]] = matrix(inits.defaults[[elem]],tmp.dim,1) }
  }
  inits$U = rbind(inits$C,inits$U) #yes, C on top
  inits$A = rbind(inits$D,inits$A) #yes, D on top
return(inits)
}