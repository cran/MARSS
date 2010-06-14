MARSSkemcheck = function(modelObj){
#This checks that the model can be handled by the MARSSkem algorithm
fixed = modelObj$fixed
free = modelObj$free
m=dim(fixed$Q)[1]; n=dim(fixed$R)[1]
TT=dim(modelObj$data)[2]
errmsg = "Try using foo=MARSS(..., fit=FALSE), then print(foo$model) to see what model you are trying to fit.\n"

  # If T=2 then kf will break
  if(TT<=2)
    stop("MARSSkemcheck: the number of time steps is <=2.  More than 2 data points are needed to estimate parameters.\n", call.=FALSE)

  # Check V0 setting and set to fixed value
  if(!is.fixed(fixed$V0))  # No part of V0 is ever estimated in our EM algorithm; either it's a prior (non-zero value) or 0 
      stop("MARSSkemcheck: fixed$V0 has NAs.  No part of V0 is ever estimated in our EM algorithm.  Use fixed$V0 to set a prior.\n", errmsg, call.=FALSE)
  if(!is.fixed(fixed$x0) &  !(identical(unname(fixed$V0), array(0,dim=c(m,m))))) 
      stop("MARSSkemcheck: part of x0 is being estimated.  In this case, fixed$V0 must equal zero so x0 is treated as fixed but unknown. Use fixed x0 and V0 to set a prior.\n", errmsg, call.=FALSE)

  constr.type = describe.marssm(modelObj)

  ############ MARSS 1.0 does not allow B to be estimated and it must be diagonal (even if fixed)
   if(!is.fixed(fixed$B) )
      stop("MARSSkemcheck: MARSS 1.0 does not allow B to be estimated.\n", call.=FALSE )
   if(!all(abs(eigen(fixed$B)$values)<=1)) 
      stop("MARSSkemcheck: in MARSS 1.0 all the eigenvalues of B must be within the unit circle: all(abs(eigen(fixed$B)$values)<=1)\n", errmsg, call.=FALSE)
              
  ############ Check that Z form does not conflict with R form
   if(!is.fixed(fixed$Z) )
      stop("MARSSkemcheck: MARSS 1.0 does not allow Z to be estimated.)\n", errmsg, call.=FALSE )

  ############ Check that R is diagonal if there are missing values
  if(any(modelObj$data==modelObj$miss.value) && !is.fixed(fixed$R) && n!=1) {
    #Then it must be diagonal
    R.first.word=substr(constr.type$R,1,8)
    if(R.first.word != "diagonal") stop("MARSSkemcheck: If there are missing values, R must be scalar, diagonal or fixed.\n", errmsg, call.=FALSE)
    }
    
return(constr.type)
}
