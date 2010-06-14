###################################################################################
# as.marssm(popWrap) function
# Coerce popWrap wrapper object to MARSS model object (class marssm).
###################################################################################
as.marssm <- function(wrapperObj)
{
  ## Check object passed in
  if (class(wrapperObj) != "popWrap") 
    stop("Argument to as.marssm() must be a popWrap object.")

  if(is.null(wrapperObj$fixed)) fixed = list()
  else fixed = wrapperObj$fixed
  if(is.null(wrapperObj$free)) free = list()
  else free = wrapperObj$free
  constraint = wrapperObj$constraint

  ## Z dealt with in popWrap() 
  # KW
  if(is.vector(wrapperObj$data)) data = matrix(wrapperObj$data, nrow=1)
  else data = wrapperObj$data 
  n = dim(data)[1]
  m = wrapperObj$m

  model.dims = list(Z=c(n,m),U=c(m,1),A=c(n,1),B=c(m,m),Q=c(m,m),R=c(n,n),x0=c(m,1),V0=c(m,m))

  ## Translate the constraints names (shortcuts) into fixed and free
  ## fixed NAs for those values that are not fixed
  ## free NAs for those values that are not free and shared numbers for those values that are shared 

    if(is.factor(constraint$Z)) { 
      X.names = unique(constraint$Z)
      fixed$Z <- matrix(0,model.dims$Z[1], model.dims$Z[2])  
      for(i in X.names) fixed$Z[which(constraint$Z==i), which(as.vector(X.names)==i)] <- 1
      free$Z = array(NA, dim=model.dims$Z)
    }
    
  for(el in model.elem.w.V0) {
    if(is.factor(constraint[[el]]) && model.dims[[el]][1]==model.dims[[el]][2] && el!="Z") {
      free[[el]]=array(NA,dim=model.dims[[el]])
      diag(free[[el]])=as.character(constraint[[el]])
      fixed[[el]]=array(0,dim=model.dims[[el]])
      diag(fixed[[el]])=NA 
    }
    if(is.factor(constraint[[el]]) && model.dims[[el]][2]==1 && el!="Z") {
      free[[el]]=array(constraint[[el]],dim=model.dims[[el]])
      fixed[[el]] = array(NA,dim=model.dims[[el]])
    }
    if( identical(constraint[[el]],"identity") ) { 
      free[[el]]=array(NA,dim=model.dims[[el]])
      fixed[[el]]=makediag(1,nrow=model.dims[[el]][1])
    } #m=n
  if(identical(constraint[[el]],"diagonal and equal")) {
    free[[el]] = array(NA, dim=model.dims[[el]])
    diag(free[[el]]) = 1
    fixed[[el]] = array(0, dim=model.dims[[el]])
    diag(fixed[[el]]) = NA
  }
  if(identical(constraint[[el]],"diagonal and unequal")) {
    free[[el]]=array(NA,dim=model.dims[[el]])
    diag(free[[el]])=1:model.dims[[el]][1]
    fixed[[el]]=array(0,dim=model.dims[[el]])
    diag(fixed[[el]])=NA 
  }
  if(identical(constraint[[el]],"unconstrained") || identical(constraint[[el]],"unequal")) {
    free[[el]]=array(seq(1,model.dims[[el]][1]*model.dims[[el]][2]),dim=model.dims[[el]])
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  }
  if(identical(constraint[[el]],"equalvarcov")) {
    free[[el]]=makediag(rep(1,model.dims[[el]][1]),nrow=model.dims[[el]][1])+array(1,dim=model.dims[[el]])
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  }
    if(identical(constraint[[el]],"equal")) { 
    free[[el]]=array(1,dim=model.dims[[el]])
    fixed[[el]] = array(NA,dim=model.dims[[el]]) 
  }
  if(identical(constraint[[el]],"zero")) { 
    free[[el]]=array(NA,dim=model.dims[[el]])
    fixed[[el]] = array(0,dim=model.dims[[el]]) 
  }
  if(identical(constraint[[el]],"ones")) { 
    free[[el]]=array(NA,dim=model.dims[[el]])
    fixed[[el]] = array(1,dim=model.dims[[el]])
    if(el %in% c("Q","R","B")) fixed[[el]] = makediag(1,nrow=model.dims[[el]][1])
  }
  if(is.matrix(constraint[[el]])) {
      fixed[[el]] = constraint[[el]]
      free[[el]] = array(1:length(fixed[[el]]),dim=dim(fixed[[el]]))
      free[[el]][!is.na(fixed[[el]])] = NA
      }

  } # end for(el in model.elem)

  if(identical(constraint$A, "scaling")) {  
    ## Construct A from fixed Z matrix   
    fixed$A = matrix(NA,model.dims$A[1],model.dims$A[2])
    free$A = matrix(1:model.dims$A[1],model.dims$A[1],1)
    for(i in 1:model.dims$Z[2]) {
      fixed$A[min(which(fixed$Z[,i]==1)), 1] <- 0
      free$A[min(which(fixed$Z[,i]==1)), 1] <- NA
    }
  }   
    
  ## Create free matrices when not passed in (b509)
  el = model.elem
  nofree = el[!(el %in% names(free)) & (el %in% names(fixed))]
  for (i in nofree) { #any NAs will be estimated separately
      free[[i]] = array(1:length(fixed[[i]]),dim=dim(fixed[[i]]))
      free[[i]][!is.na(fixed[[i]])] = NA
    }
    
  ## Create fixed matrices when not passed in and no NAs in free (b509)
  el = model.elem
  nofixed = el[!(el %in% names(fixed)) & (el %in% names(free))]
  for (i in nofixed) {
    if(!any(is.na(free[[i]]))) fixed[[i]] = array(NA, dim=dim(free[[i]])) 
    } 

  ## Make free character to avoid problems with unique(), table(), etc.
  for (i in 1:length(free)) mode(free[[i]]) = "character"

  ## Create marssm obj

  modelObj <- marssm(fixed=fixed, free=free, data=data, miss.value=wrapperObj$miss.value)

  ## Return marssm obj
  modelObj
} 


