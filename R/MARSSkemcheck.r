MARSSkemcheck = function(modelObj, method="kem", kf.x0=NULL){
#This checks that the model can be handled by the MARSSkem algorithm
fixed = modelObj$fixed
free = modelObj$free
m=dim(fixed$Q)[1]; n=dim(fixed$R)[1]
TT=dim(modelObj$data)[2]
errmsg = " Try using foo=MARSS(..., fit=FALSE), then summary(foo$model) to see what model you are trying to fit.\n"

#ensure that kf.x0 is passed in
  if(is.null(kf.x0))
    stop("Stopped in MARSSkemcheck() because kf.x0 not passed into function.\n", call.=FALSE)

# If T=2 then kf will break
  if(TT<=2)
    stop("Stopped in MARSSkemcheck() because the number of time steps is <=2.\nMore than 2 data points are needed to estimate parameters.\n", call.=FALSE)

  ############ Check that B is within the unit circle
   if(is.fixed(fixed$B) && !all(abs(eigen(fixed$B,only.values=TRUE)$values)<=1)){ 
      msg=c(" In MARSS 2.0 all the eigenvalues of B must be within the unit circle: all(abs(eigen(fixed$B)$values)<=1)\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }

   ############ Check that if R has 0s, then the corresponding row of A and Z are fixed
   # and that if x0 is fixed, data = Z%*%x0
   diag.R=takediag(fixed$R);  diag.R[is.na(diag.R)]=1
   R.degen.elements=(diag.R==0)
   if(any(R.degen.elements)){
      allowed.to.be.degen = !is.na(fixed$A[R.degen.elements]) & !apply(is.na(fixed$Z[R.degen.elements,,drop=FALSE]),1,any)
      if(any(!allowed.to.be.degen,na.rm=TRUE)) {
        msg=c(" If an element of the diagonal of R is 0, the corresponding row of A and Z must be fixed.\n", errmsg)
        cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="")
        stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }
      
      if(kf.x0=="x10"){
        fixed.hat.y0=fixed$Z[R.degen.elements,,drop=FALSE]%*%fixed$x0+fixed$A[R.degen.elements,drop=FALSE]
        data1 = modelObj$data[R.degen.elements,1]
        data1[data1==modelObj$miss.value]=NA #replace miss.value with NA
        both.fixed=!is.na(fixed.hat.y0) & !is.na(data1)
        if(any(both.fixed) & any(!(data1[both.fixed]==fixed.hat.y0[both.fixed]))){
          msg=c(" If an element of the diagonal of R is 0, the corresponding A, Z, and x0 are fixed, \n
            and kf.x0=x10 and data are not missing at t=1, then the data must match the fixed values.\n", errmsg)
          cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="")
          stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
        }
      }
   }
   
   ############ Check that R and Q don't have 0s in the same place 
   #Z where 0=0 and !0=1
   Z.01 = fixed$Z; Z.01[is.na(Z.01)]=1; Z.01=(Z.01!=0)
   diag.Q=diag(fixed$Q);  diag.Q[is.na(diag.Q)]=TRUE; diag.Q[diag.Q!=0]=TRUE; diag.Q[diag.Q==0]=FALSE
   Q.degen.elements=!(Z.01%*%diag.Q)
   if(any(Q.degen.elements & R.degen.elements)) {
      msg=c(" Warning: an element of the diagonal of R is 0, and the corresponding row of Q is also 0.\n This might lead to MARSS errors.\n", errmsg)
      cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
      #stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      }
  
  ############ Check that the B sub matrix for Q=0 is diagonal  
    diag.Q=takediag(modelObj$fixed$Q)
    if( any(diag.Q==0,na.rm=TRUE) ){
       B.0 = fixed$B[diag.Q==0, diag.Q==0]
       if(!is.diagonal(B.0, na.rm=TRUE)){ 
          msg=c(" The B sub matrix, corresponding to Q diagonal = 0, must be diagonal.\n", errmsg)
          cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
          stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
       }
    }    

constr.type = describe.marssm(modelObj)   
return(constr.type)
}
