#######################################################################################################
#   MARSSkf function
#   Utility function to choose the Kalman filter and smoother
#######################################################################################################
MARSSkf = function( MLEobj, only.logLik=FALSE, return.lag.one=TRUE, return.kfas.model=FALSE ) {
if(MLEobj$fun.kf=="MARSSkfss") 
  return(MARSSkfss(MLEobj))
if(MLEobj$fun.kf=="MARSSkfas") 
  return( MARSSkfas(MLEobj, only.logLik=FALSE, return.lag.one=TRUE, return.kfas.model=FALSE ) )
return(list(ok=FALSE, errors="kf.function does not specify a valid Kalman filter and smoother function."))
}