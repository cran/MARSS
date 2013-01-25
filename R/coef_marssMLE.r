###############################################################################################################################################
#  coef method for class marssMLE. 
##############################################################################################################################################
coef.marssMLE <- function (object, ..., type="list", form=NULL) {
#First make sure specified equation form has a corresponding function to do the conversion to marssm object
  return.obj=list()
  what.to.print = type
  if(is.null(form)){ #allow user to printing using a different form than call$form
    if(is.null(object[["call"]][["form"]])) form="marss" else form=object$call$form
  }
  coef.fun = paste("coef_",form,sep="")
  tmp=try(exists(coef.fun,mode="function"),silent=TRUE)
  if(isTRUE(tmp)){
      #the coef function returns an updated object based on form
      object=eval(call(coef.fun, object))
    }
  if(class(object)=="marssMLE"){ #use object

   fixed=object$model$fixed
    free=object$model$free

if(!(type %in% c("vector", "list", "par", "matrix", names(fixed))))
  stop("coef.marssMLE: allowed coef types are vector, list, par, matrix")
  
  for(what in what.to.print){    
    if(what=="par" | what=="list"){
      return.obj[[what]]=object$par
    }
    if(what=="vector"){
      paramvector = NULL
      for(elem in names(fixed)){
        if(dim(object$par[[elem]])[1]>0){ #there are estimates
        mat.names = colnames(free[[elem]])
        tmp = as.vector(object$par[[elem]]) 
        mat.names = paste(rep(elem, length(mat.names)), rep(".", length(mat.names)), mat.names, sep="")
        names(tmp) = mat.names
        paramvector = c(paramvector, tmp)
        }
      }
      return.obj[[what]]=paramvector
    }

    #if there is info on the matrix dimensions in the model list, use it
    if(is.null(object[["model"]][["model.dims"]])) stop("In coef.marssMLE: model list is missing model.dims")
    par.dims=object$model$model.dims
    if(what == "matrices" | what=="matrix"){
      par.mat=list()
      for(elem in names(fixed)){
        par.dim = par.dims[[elem]]
        if(length(par.dim)!=2) stop("In coef.marssMLE: par.dim is not a length 2 vector")
        the.par = parmat(object,elem,t=1:max(dim(fixed[[elem]])[3],dim(free[[elem]])[3]), dims=par.dim)[[elem]]
        par.mat[[elem]]=the.par
      }
      return.obj[[what]]=par.mat
    }    
    if(what %in% names(fixed)){
      the.par = parmat(object,what,t=1:max(dim(fixed[[what]])[3],dim(free[[what]])[3]), dims=par.dim)[[what]]
      return.obj[[what]]=the.par
    }
    if(what %in% names(object$model)){
      return.obj[[what]]=object$model[[what]]
    }
  } #for what in what.to.print
  if(length(return.obj)==0) return.obj=NULL
  if(length(return.obj)==1) return.obj=return.obj[[1]]
  return(return.obj)
  }else{ #object is  not type marssMLE
  stop("coef.marssMLE: coef needs a marssMLE object.")
  }   
 }  #end of coef.marssMLE
 
 coef_marss = function(x){ return(x) }