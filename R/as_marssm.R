###################################################################################
# as.marssm(popWrap) function
# Coerce popWrap wrapper object to MARSS model object (class marssm).
###################################################################################
as.marssm <- function(wrapperObj)
{
  ## Check object passed in
  if (class(wrapperObj) != "popWrap") 
    stop("Stopped in as.marssm() because argument to as.marssm() must be a popWrap object.\n", call.=FALSE)

  if(is.null(wrapperObj$fixed)) fixed = list()
  else fixed = wrapperObj$fixed
  if(is.null(wrapperObj$free)) free = list()
  else free = wrapperObj$free
  model = wrapperObj$model

  ## Z dealt with in popWrap() 
  # KW
  if(is.vector(wrapperObj$data)) data = matrix(wrapperObj$data, nrow=1)
  else data = wrapperObj$data 
  n = dim(data)[1]
  m = wrapperObj$m

  model.dims = list(Z=c(n,m),U=c(m,1),A=c(n,1),B=c(m,m),Q=c(m,m),R=c(n,n),x0=c(m,1),V0=c(m,m))

  ## Translate the model structure names (shortcuts) into fixed and free
  ## fixed NAs for those values that are not fixed
  ## free NAs for those values that are not free and shared numbers for those values that are shared 

    if(is.factor(model$Z)) { 
      X.names = unique(model$Z)
      fixed$Z <- matrix(0,model.dims$Z[1], model.dims$Z[2])  
      for(i in X.names) fixed$Z[which(model$Z==i), which(as.vector(X.names)==i)] <- 1
      free$Z = array(NA, dim=model.dims$Z)
    }
    if(identical(model$Z,"onestate")) { 
      free$Z=array(NA,dim=model.dims[[el]])
      fixed$Z = array(1,dim=model.dims[[el]])
    }
  for(el in model.elem) {
    if(is.factor(model[[el]]) && model.dims[[el]][1]==model.dims[[el]][2] && el!="Z") {
      free[[el]]=array(NA,dim=model.dims[[el]])
      diag(free[[el]])=as.character(model[[el]])
      fixed[[el]]=array(0,dim=model.dims[[el]])
      diag(fixed[[el]])=NA 
    }
    if(is.factor(model[[el]]) && model.dims[[el]][2]==1 && el!="Z") {
      free[[el]]=array(model[[el]],dim=model.dims[[el]])
      fixed[[el]] = array(NA,dim=model.dims[[el]])
    }
    if( identical(model[[el]],"identity") ) { 
      free[[el]]=array(NA,dim=model.dims[[el]])
      fixed[[el]]=makediag(1,nrow=model.dims[[el]][1])
    } #m=n
  if(identical(model[[el]],"diagonal and equal")) {
    free[[el]] = array(NA, dim=model.dims[[el]])
    diag(free[[el]]) = 1
    fixed[[el]] = array(0, dim=model.dims[[el]])
    diag(fixed[[el]]) = NA
  }
  if(identical(model[[el]],"diagonal and unequal")) {
    free[[el]]=array(NA,dim=model.dims[[el]])
    diag(free[[el]])=1:model.dims[[el]][1]
    fixed[[el]]=array(0,dim=model.dims[[el]])
    diag(fixed[[el]])=NA 
  }
  if(identical(model[[el]],"unconstrained")){
  if(el %in% c("Q","R","V0")){  #variance-covariance matrices
    dim.mat = model.dims[[el]][1]
    free[[el]]=array(NA,dim=model.dims[[el]]) 
    for(i in 1:dim.mat){
      free[[el]][i,i]=paste("(",i,",",i,")",sep="")
      for(j in 1:dim.mat) free[[el]][i,j]=free[[el]][j,i]=paste("(",i,",",j,")",sep="")
    }
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  }else{ #not var-cov matrix
    free[[el]]=array(seq(1,model.dims[[el]][1]*model.dims[[el]][2]),dim=model.dims[[el]])
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  } }
  if( identical(model[[el]],"unequal")) {
    free[[el]]=array(seq(1,model.dims[[el]][1]*model.dims[[el]][2]),dim=model.dims[[el]])
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  }
  if(identical(model[[el]],"equalvarcov")) {
    free[[el]]=makediag(rep(1,model.dims[[el]][1]),nrow=model.dims[[el]][1])+array(1,dim=model.dims[[el]])
    fixed[[el]]=array(NA,dim=model.dims[[el]])
  }
    if(identical(model[[el]],"equal")) { 
    free[[el]]=array(1,dim=model.dims[[el]])
    fixed[[el]] = array(NA,dim=model.dims[[el]]) 
  }
  if(identical(model[[el]],"zero")) { 
    free[[el]]=array(NA,dim=model.dims[[el]])
    fixed[[el]] = array(0,dim=model.dims[[el]]) 
  }
  if(is.matrix(model[[el]])) {
      if(is.numeric(model[[el]])){
        fixed[[el]] = model[[el]]
        free[[el]] = array(NA,dim=dim(fixed[[el]]))
        if(sum(is.na(fixed[[el]]))!=0) free[[el]][is.na(fixed[[el]])] = 1:sum(is.na(fixed[[el]]))
      }
      if(is.character(model[[el]])){
        free[[el]] = model[[el]]
        fixed[[el]] = array(NA,dim=dim(free[[el]]))
      }
      if(is.list(model[[el]])){
        mat = model[[el]]
        dim.mat = dim(mat)
        free[[el]] = array(NA,dim=dim.mat)
        free[[el]][array(sapply(mat,is.character),dim=dim.mat)]=unlist(mat[array(sapply(mat,is.character),dim=dim.mat)])
        fixed[[el]] = array(NA,dim=dim.mat)
        fixed[[el]][array(sapply(mat,is.numeric),dim=dim.mat)]=unlist(mat[array(sapply(mat,is.numeric),dim=dim.mat)])
      }
      }

  } # end for(el in model.elem)

  #Set the column names on Z; otherwise this information is lost
    if(is.matrix(model$Z)) colnames(fixed$Z)=colnames(model$Z)
    if(is.factor(model$Z)) colnames(fixed$Z)=unique(model$Z)

  if(identical(model$A, "scaling")) {  
    ## Construct A from fixed Z matrix   
    fixed$A = matrix(NA,model.dims$A[1],model.dims$A[2])
    free$A = matrix(1:model.dims$A[1],model.dims$A[1],1)
    for(i in 1:model.dims$Z[2]) {
      fixed$A[min(which(fixed$Z[,i]==1)), 1] <- 0
      free$A[min(which(fixed$Z[,i]==1)), 1] <- NA
    }
  }   
    
  ## Create free matrices when not passed in (b509) but fixed is passed in
  el = model.elem
  nofree = el[!(el %in% names(free)) & (el %in% names(fixed))]
  for (i in nofree) { #any NAs will be estimated separately
      free[[i]] = array(NA,dim=dim(fixed[[i]]))
      if(sum(is.na(fixed[[i]]))!=0) free[[i]][is.na(fixed[[i]])] = 1:sum(is.na(fixed[[i]]))
    }
    
  ## Create fixed matrices when not passed in and no NAs in free (b509)
  el = model.elem
  nofixed = el[!(el %in% names(fixed)) & (el %in% names(free))]
  for (i in nofixed) {
    if(!any(is.na(free[[i]]))) fixed[[i]] = array(NA, dim=dim(free[[i]])) 
    } 

  ## Make free character to avoid problems with unique(), table(), etc.
  for (i in 1:length(free)) mode(free[[i]]) = "character"

  ## Ensure that R doesn't interpret NA as logical
  if(is.na(wrapperObj$miss.value)) wrapperObj$miss.value = as.numeric(wrapperObj$miss.value)

  ## Create marssm obj

  modelObj <- marssm(fixed=fixed, free=free, data=data, miss.value=wrapperObj$miss.value)

  ## Return marssm obj
  modelObj
} 


