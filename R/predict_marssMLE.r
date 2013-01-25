###############################################################################################################################################
#  Predict method for class marssMLE. 
##############################################################################################################################################
predict.marssMLE <- function (x, n.ahead=1, t.start=NULL, newdata=list(), se.fit=TRUE, form=NULL, ...) {
#newdata must reflect the inputs allowed for a specific form; y (data) is an input.
#n.ahead = 0 means to predict values WITHIN 1:TT.  n.ahead > 0 means predict (TT+1):(TT+n.ahead

  #if n.ahead is not passed in AND newdata is passed in, n.ahead will be taken from the info in newdata
  if(missing(n.ahead) & length(newdata)!=0) n.ahead=NULL
  #By default, start the predict at the end of the time series
  if(missing(t.start)) dim(x$model$data)[2]
  
  #First make sure specified equation form has a corresponding function to do the conversion to marssm object
  return.obj=list()
  orig.x=x
  if(is.null(form)){ #allow user to predict using a different form than x$call$form
    if(is.null(x[["call"]][["form"]])) form="marss" else form=x$call$form
  }
  predict.fun = paste("predict_",form,sep="")
  tmp=try(exists(predict.fun,mode="function"),silent=TRUE)
  if(isTRUE(tmp)){
      #the predict function returns an updated x with the model object corresponding to the form
      x=eval(call(predict.fun, x, newdata))
    }
  if(class(x)=="marssMLE"){
    modelObj=x$model
    ## Check that the marssm object is ok
    ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
    tmp = is.marssm(modelObj)
    if(!isTRUE(tmp)) {
      return(x) 
      stop("Stopped in predict() due to problem(s) with model specification. \nThe marssMLE object constructed for predict call is being returned.  Examine $model for problems.\n", call.=FALSE)
    }    
  }else{ #object is  not type marssMLE
    stop("predict.marssMLE: predict needs a marssMLE object.")
  } #x is type marssMLE  
 }  #end of print.marssMLE

predict_marssm = function(x){ return(x) }