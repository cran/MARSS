########################################################################
# is.marssm function
# Check that the modelObj object has all the parts it needs
# data, fixed, free, miss.value, and X.names
# and that these have the proper size and form
# m is pulled from fixed$x0
########################################################################
is.marssm <- function(modelObj) 
{
  if(class(modelObj) != "marssm") stop("Stopped in is.marssm() because object class is not marssm.\n", call.=FALSE)
  msg = NULL

###########################
# First do some basic mode and presence tests so that the other tests work
###########################

  ## Check for required components
  el = c("data","fixed","free","miss.value","X.names","tinitx","diffuse")
  if( !all(el %in% names(modelObj)) ) { 
    msg = c(msg, "Element", el[!(el %in% names(modelObj))], "is missing from the model object.\n")
  }
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssm()\n", msg)
    return(msg)
  }
  
  ## Check that free and fixed are numeric matrices with no NA or Infs
  for(mat in c("fixed","free")){
      if (!is.list(modelObj[[mat]])) msg = c(msg, paste("modelObj$",mat," must be a list of matrices.\n",sep="")) 
      for (i in 1:length(modelObj[[mat]])) {
        if(class(modelObj[[mat]][[i]]) != "array" || length(dim(modelObj[[mat]][[i]]))!=3){ 
            msg = c(msg, paste("modelObj$",mat,"$",names(modelObj[[mat]])[i]," must be a 3D matrix.\n", sep=""))
        }
        if(mode(modelObj[[mat]][[i]]) != "numeric" || any(is.na(modelObj[[mat]][[i]])) || any(is.infinite(modelObj[[mat]][[i]])) ) 
            msg = c(msg, paste("modelObj$",mat,"$",names(modelObj[[mat]])[i]," must be numeric, have no NAs, and no Infs.\n", sep=""))  
      }
  }
  
  if( length(dim(modelObj$data)) != 2)
     msg = c(msg, "Data is not a 2D matrix.\n")
       ## check for T=1
  if( !is.numeric(modelObj$data ) ) msg = c(msg, "Data must be numeric.\n")
  if( dim(modelObj$data)[2] == 1 ) msg = c(msg, "Data has only one time point.\n")

  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssm()\n", msg)
    return(msg)
  }

###########################
# Check that fixed and free are complete and consistent
###########################

  ## Set up element names that appear as fixed/free pair
  en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")
  ## Check presence of all fixed/free pairs and that they have dim=3 
  fixed=modelObj$fixed
  free=modelObj$free
  fix = est =  NULL
  for (elem in en) {
    fixed.null.flag <- is.null(fixed[[elem]]) 
    free.null.flag <- is.null(free[[elem]])
    fix <- c(fix, fixed.null.flag)
    est <- c(est, free.null.flag)
  }  
  if (any(c(fix, est))) {  #stop now since the rest of the tests won't work
    if(any(fix)) {
      msg = c(msg, paste("Missing fixed", en[fix],"\n"))   }
    if(any(est)) {
      msg = c(msg, paste("Missing free", en[est],"\n"))    }
    msg=c("\nErrors were caught in is.marssm()\n", msg)
    return(msg)
  }
 
  dim.fixed = dim.free =  NULL
  for (elem in en) {
    dim.fixed.flag = dim.free.flag = FALSE
    if(length(dim(free[[elem]]))!=3){ dim.free.flag = TRUE } #3-dimensions
    if(length(dim(fixed[[elem]]))!=3){ dim.fixed.flag = TRUE }
    dim.fixed <- c( dim.fixed, dim.fixed.flag )
    dim.free <- c( dim.free, dim.free.flag )
  }  
  if (any(c(dim.fixed, dim.free))) {  #stop now since the rest of the tests won't work
    if(any(dim.fixed)) {
      msg = c(msg, paste("fixed", en[fix],"is not 3D.\n"))   }
    if(any(dim.free)) {
      msg = c(msg, paste("free", en[est],"is not 3D.\n"))    }
    msg=c("\nErrors were caught in is.marssm()\n", msg)
    return(msg)
  }

  #Now we can check the dimensions of the fixed and free pairs
  n = dim(modelObj$data)[1]
  TT = dim(modelObj$data)[2]
  m = dim(modelObj$fixed$x0)[1]
  correct.dim1 = c(Z=n,A=n,R=n,B=m, U=m, Q=m, x0=m, V0=m)
  correct.dim2 = c(Z=m,A=1,R=n,B=m, U=1, Q=m, x0=1, V0=m)

  fix = est = dim.fixed = dim.free =  NULL

  for (elem in en) {
  ## Check for problems in the fixed/free pairs. Problems show up as TRUE 
     dim.fixed.flag = dim.free.flag = nomatch.flag = FALSE
        
     # check dim
     dim.fixed.flag = !isTRUE(all.equal( dim(fixed[[elem]])[1], correct.dim1[[elem]]*correct.dim2[[elem]] ) )
     dim.fixed.flag = dim.fixed.flag | !isTRUE( all.equal( dim(fixed[[elem]])[2], 1 ) )
     #test that dim3 is either 1 or TT
     dim.fixed.flag = dim.fixed.flag | (!isTRUE( all.equal( dim(fixed[[elem]])[3], TT ) ) & !isTRUE(all.equal( dim(fixed[[elem]])[3], 1 ) )) 
     dim.free.flag = !isTRUE(all.equal( dim(free[[elem]])[1], correct.dim1[[elem]]*correct.dim2[[elem]] ) )
     dim.free.flag = dim.free.flag | (!isTRUE(all.equal( dim(free[[elem]])[3], TT ) ) & !isTRUE(all.equal( dim(free[[elem]])[3], 1 ) ))

    dim.fixed = c(dim.fixed, dim.fixed.flag)
    dim.free = c(dim.free, dim.free.flag) 
   }  
  if (any(c(dim.fixed, dim.free))) {  #There's a problem
    if(any(dim.fixed)) {
      msg = c(msg, paste("fixed", en[dim.fixed], "dims are incorrect. Dims should be (", correct.dim1[dim.fixed], "x", correct.dim2[dim.fixed],", 1) based on data and fixed$x0.\n"))
    }
    if(any(dim.free)) {
      msg = c(msg, paste("free", en[dim.free], "dims are incorrect. Dim 1 be ", correct.dim1[dim.free], "x", correct.dim2[dim.free],"based on data and fixed$x0.\n"))
    }
   msg=c("\nErrors were caught in is.marssm()\n", msg)
   return(msg)
  } 

###########################
# Check that x0 and V0 are not time-varying
###########################
  en = c("x0", "V0")
  time.var = NULL
  for (elem in en) {
    time.var.flag = FALSE
    time.var.flag = dim(fixed[[elem]])[3]!=1 | dim(free[[elem]])[3]!=1
    time.var <- c(time.var, time.var.flag)
  }
  if(any(time.var)) {  #There's a problem
    msg = c(msg, paste(en[time.var], "cannot be time-varying.  3rd dim of fixed and free must equal 1.\n"))
    msg=c("\nErrors were caught in is.marssm()\n", msg)
    return(msg)
  }
   
###########################
# Check that none of the var-cov matrices have negative values on the diagonal
# and that there are no f+Dq elements only f+0q or 0+Dq
# and D must be a design matrix, so no beta_1*q1 + beta_2*q2 elements
###########################
  en = c("R", "Q", "V0")
  neg = bad.var = not.design = NULL
  for (elem in en) {
    neg.flag = bad.var.flag = not.design.flag = FALSE
    for(i in 1:max(dim(free[[elem]])[3],dim(fixed[[elem]])[3])){
      if(dim(fixed[[elem]])[3]==1){i1=1}else{i1=i}
      if(dim(free[[elem]])[3]==1){i2=1}else{i2=i}
      if(is.fixed(free[[elem]][,,min(i,dim(free[[elem]])[3]),drop=FALSE])){ #this works on 3D mats
        zero.free.rows = matrix(TRUE,correct.dim1[[elem]]*correct.dim2[[elem]],1)
      }else{
        zero.free.rows=apply(free[[elem]][,,i2,drop=FALSE]==0,1,all) #works on 3D mat
        #the requirement is that each estimated element (in p) appears only in one place in the varcov mat, but fixed rows (0 rows) are ok
        not.design.flag = !is.design(free[[elem]][,,i2,drop=FALSE],strict=FALSE,zero.rows.ok=TRUE,zero.cols.ok=TRUE) #works on 3D if dim3=1
      }
      zero.fixed.rows=apply(fixed[[elem]][,,i1,drop=FALSE]==0,1,all) #works on 3D
      fixed.mat = unvec(fixed[[elem]][,,i1],dim=c(correct.dim1[[elem]],correct.dim2[[elem]]))
      if( any(!zero.fixed.rows & !zero.free.rows) ) bad.var.flag = TRUE   #no f+Dq rows
      if(any(takediag(fixed.mat)<0,na.rm=TRUE)) neg.flag=TRUE      #no negative diagonals
    } #end the for loop over time
    not.design = c(not.design, not.design.flag)
    neg = c(neg, neg.flag)
    bad.var = c(bad.var, bad.var.flag)
    } #enf the for loop over elem
  if(any(neg)) {
      msg = c(msg, paste("Negative values are on the diagonal of ", en[neg], ". Neg values are illegal on the diag of a var-cov matrix.\n", sep=""))
      }
  if(any(bad.var)) {
      msg = c(msg, paste("Fixed and estimated values are combined in some elements of ", en[bad.var], ". This is not allowed.\n", sep=""))
      }
  if(any(not.design)) {
      msg = c(msg, paste("The D matrices of ", en[not.design], " must be design matrices.\n", sep=""))
      }
      
###########################
# TURNED OFF Check that the Q and R matrices have no zeros
###########################
if(1==0){ #TURNED OFF
  en = c("R", "Q")
  #this one catches if the user set 1 0 on the diagonal and set the others to NA (estimated)
  zer = NULL
  for (elem in en) {
    zero.flag = FALSE
    if(any(takediag(modelObj$fixed[[elem]])==0,na.rm=TRUE)) zero.flag=TRUE
    zer = c(zer, zero.flag)
    }
  if(any(zer)) { 
      msg = c(msg, paste("Zeros are on the diagonal of ", en[zer], ". Zeros are illegal on the diag of the R or Q matrices.\n", sep=""))
      }
}

###########################
# Check that fixed V0, Q and R matrices are symmetric and positive-definite
###########################
  en = c("R", "Q", "V0")
  pos = symm = NULL
  for (elem in en) {
    symm.flag = FALSE
    pos.flag = FALSE
    var.dim = c(correct.dim1[[elem]],correct.dim2[[elem]])
    for(i in 1:max(dim(fixed[[elem]])[3],dim(free[[elem]])[3])){
      if(dim(fixed[[elem]])[3]==1){i1=1}else{i1=i}
      if(dim(free[[elem]])[3]==1){i2=1}else{i2=i}
      #works on 3D if dim3=1
      par.as.list = fixed.free.to.formula(fixed[[elem]][,,i1,drop=FALSE],free[[elem]][,,i2,drop=FALSE],var.dim) #coverts the fixed,free pair to a list matrix
      if(!isTRUE(all.equal(par.as.list, t(par.as.list)))) symm.flag=TRUE
      if(is.fixed(free[[elem]][,,i2,drop=FALSE])){
        var.mat=unvec(fixed[[elem]][,,i1,drop=FALSE],dim=var.dim)
        tmp = try( eigen(var.mat), silent=TRUE )
        if(class(tmp)=="try-error") pos.flag=TRUE
        else if(!all(tmp$values >= 0)) pos.flag=TRUE
      }
    } #end for loop over time
    pos = c(pos, pos.flag)
    symm = c(symm, symm.flag)
  } #end for loop over elements
  if(any(pos)) msg = c(msg, paste("The fixed matrix ", en[pos], " is not positive-definite (and var-cov matrices must be).\n", sep=""))
  if(any(symm)) msg = c(msg, paste("The variance matrix ", en[symm], " is not symmetric (and var-cov matrices must be).\n", sep=""))
        
###########################
# Check data and missing values consistency if data present
###########################
    if(!is.numeric(modelObj$data)) msg = paste(msg, "Data must be numeric. \n")
    for( bad.val in c(NA, NaN, Inf, -Inf)){
      if(!identical(bad.val, modelObj$miss.value) && ( bad.val %in% modelObj$data ) ){  
            msg = c(msg, paste("Data cannot have ", bad.val, "s unless miss.value is ", bad.val,". \n",sep="")) }
      }

      ## Check that data are not all missing values 
      if(is.na(modelObj$miss.value)){
        if(all(is.na(modelObj$data))){ msg = c(msg, "Data consists only of missing values.\n") }
        }
        else {
          if(all(modelObj$data == modelObj$miss.value) ){ msg = c(msg, "Data consists only of missing values.\n") }
        }

###########################
# Check X.names
###########################
    if(length(modelObj$X.names) != m){ msg = c(msg, "Length of X.names does not match model.\n") }
  
###########################
# Check tinitx; must be 0 or 1
###########################
    if( !(modelObj$tinitx %in% c(0,1)) ){ msg = c(msg, "tinitx (t associated with initial x) must be 0 or 1.\n") }

###########################
# Check diffuse; must be TRUE or FALSE
###########################
    if( !(modelObj$diffuse %in% c(FALSE, TRUE)) ){ msg = c(msg, "diffuse must be TRUE or FALSE.\n") }

if(length(msg) == 0){ return(TRUE)
}else {
  msg=c("\nErrors were caught in is.marssm()\n", msg)
  return(msg)
}
}
