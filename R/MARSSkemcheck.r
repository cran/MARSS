MARSSkemcheck = function(modelObj, method="kem"){
#This checks that the model can be handled by the MARSSkem algorithm
fixed = modelObj$fixed
free = modelObj$free
m=dim(fixed$Q)[1]; n=dim(fixed$R)[1]
TT=dim(modelObj$data)[2]
errmsg = " Try using foo=MARSS(..., fit=FALSE), then print(foo$model) to see what model you are trying to fit.\n"

  # If T=2 then kf will break
  if(TT<=2)
    stop("Stopped in MARSSkemcheck() because the number of time steps is <=2.\nMore than 2 data points are needed to estimate parameters.\n", call.=FALSE)

  # Check V0 setting and set to fixed value
  if(!is.fixed(fixed$V0)){  # No part of V0 is ever estimated in our EM algorithm; either it's a prior (non-zero value) or 0 
      msg=c(" fixed$V0 has NAs.  No part of V0 is ever estimated in our EM algorithm.  Use fixed$V0 to set a prior.\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
  }
  if(!is.fixed(fixed$x0) &  !(identical(unname(fixed$V0), array(0,dim=c(m,m))))){ 
      msg=c(" Part of x0 is being estimated.  In this case, fixed$V0 must equal zero so x0 is treated as fixed but unknown. Use fixed x0 and V0 to set a prior.\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }

  constr.type = describe.marssm(modelObj)

  ############ MARSS 1.0 does not allow B to be estimated and it must be diagonal (even if fixed)
   if(!is.fixed(fixed$B) && method=="kem"){
      msg=" MARSS 1.0 does not allow B to be estimated.\n"
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }
   if(method=="kem" && !all(abs(eigen(fixed$B)$values)<=1)){ 
      msg=c(" In MARSS 1.0 all the eigenvalues of B must be within the unit circle: all(abs(eigen(fixed$B)$values)<=1)\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }
              
  ############ Check that Z form does not conflict with R form
   if(!is.fixed(fixed$Z) ){
      msg=c(" MARSS 1.0 does not allow Z to be estimated.\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }

  ############ Check that R is diagonal if there are missing values
  if((modelObj$miss.value %in% modelObj$data) && !is.fixed(fixed$R) && n!=1) {
    #Then it must be diagonal
    R.first.word=substr(constr.type$R,1,8)
    if(R.first.word != "diagonal"){
      msg=c(" If there are missing values, R must be scalar, diagonal or fixed.\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
     }
    }
    
return(constr.type)
}
