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
  
  diag.Q=takediag(modelObj$fixed$Q)
  if( any(diag.Q==0,na.rm=TRUE) ){
    B.0 = fixed$B[diag.Q==0, ,drop=FALSE]
    ############ Check that the B sub matrix for Q=0 is fixed  
    if( any( is.na(B.0) ) ){
          msg=c(" All rows of B corresponding to Q diagonal = 0 must be fixed values.\n", errmsg)
          cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
          stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
    }    

   for(elem in c("U","x0")){
   elem.0 = fixed[[elem]][diag.Q==0, ,drop=FALSE]
   ############ If u^(0) is estimated, B.0 must be independent of B.plus for the estimated u.0s 
   if( any( is.na(elem.0) ) && any(diag.Q!=0,na.rm=TRUE) ){
    B.00 = fixed$B[is.na(fixed[[elem]]) & diag.Q==0, diag.Q==0]
    B.0.plus = fixed$B[is.na(fixed[[elem]]) & diag.Q==0, diag.Q!=0]
    B.plus.0 = fixed$B[diag.Q!=0, is.na(fixed[[elem]]) & diag.Q==0]
    block.B = all( B.0.plus == 0 ) && all( B.plus.0 == 0 )
    zero.B00 =  all( B.00 == 0 ) && all(is.na(rowSums( B.0.plus )) | rowSums( B.0.plus )>0 )
    ############ Check that the abs of all eigenvalues of the B sub matrix for Q=0 are less than or = to 1, B is block diag  
    if( block.B ){
        tmp=eigen(B.00, only.values = TRUE)$values
    if( any( tmp > 1 ) ){
          msg=c(paste(" There are 0s on the diagonal of Q and corresponding ",elem,"'s are estimated.\n In this case, the abs of all eigenvalues of B.0 block must be less than or = 1.\n",sep=""), errmsg)
          cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
          stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
      } 
    }
    if( !block.B & !zero.B00 ){
          msg=c(paste(" There are 0s on the diagonal of Q and corresponding ", elem,"'s are estimated.\n In this case, B must be block diagonal (B.0 block and B.plus block)\n or B00 must be all 0s.\n",sep=""), errmsg)
          cat("\n","Errors were caught in MARSSkemcheck \n", msg, sep="") 
          stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE)
    } 
   }
   } 
} #zeros on the diagonal of Q

constr.type = describe.marssm(modelObj)   
return(constr.type)
}
