########################################################################
# marssm function
# Creates MARSS model object (class marssm).
########################################################################
marssm <- function(data=NULL, fixed, free, miss.value=-99)
{
  M <- NULL
  
  # Construct the M array for handling missing values
  # M%*%Z zeros out the row of Z when data missing 
  if (!is.null(data)) {
    n = dim(data)[1]
    TT = dim(data)[2]
    M <- array(0, dim=c(n, n, TT))  
    for(i in 1:TT){ 
      if(is.na(miss.value)){ tmp = !is.na(data[,i])
      }else { tmp= data[,i]!=miss.value }
      M[,,i] <- makediag(ifelse(tmp,1,0), nrow=n) #ifelse turns T/F into 1/0
      }
  }

  ## Make free character to avoid problems with unique(), table(), etc.
  for (i in 1:length(free)) mode(free[[i]]) = "character"

  ## Create marssm object
  modelObj <- list(fixed=fixed, free=free, data=data, M=M, miss.value=miss.value)
  class(modelObj) <- "marssm"

  modelObj

}
