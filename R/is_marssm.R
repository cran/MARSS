########################################################################
# is.marssm function
# Check that the modelObj object has all the parts it needs 
# and that these have the proper size and form
########################################################################

is.marssm <- function(modelObj) 
{
  if(class(modelObj) != "marssm") stop("is.marssm: Object class is not marssm.")
  msg = NULL

###########################
# First do some basic mode and presence tests so that the other tests work
###########################

  ## Check for required components; at a minimum free and fixed must be present
  el = c("fixed","free")
  if( !all(el %in% names(modelObj)) ) { 
    msg = c(msg, "Element", el[!(el %in% names(modelObj))], "is missing from the model object. Fixed and free work as a pair.\n")
  }

  ## Check that free is character or all NAs to avoid problems with unique(), table(), etc.
  for (i in 1:length(modelObj$free)) {
    if(mode(modelObj$free[[i]]) != "character" && !all(is.na(modelObj$free[[i]]))) 
        msg = c(msg, paste("modelObj$free$",names(modelObj$free)[i]," must be a character matrix (NAs ok).\n", sep="")) 
  }
  if (!is.list(modelObj$free)) msg = c(msg, "modelObj$free must be a list of character matrices (NAs ok).\n") 

  ## Check that fixed is numeric or all NAs
  chk.is.num = TRUE
  for (i in 1:length(modelObj$fixed)) {
    if(mode(modelObj$fixed[[i]]) != "numeric"  && !all(is.na(modelObj$fixed[[i]]))) 
            msg = c(msg, paste("modelObj$fixed$",names(modelObj$fixed)[i]," must be a numeric matrix (NAs ok).\n", sep="")) 
  }
  if (!is.list(modelObj$fixed)) msg = c(msg, "modelObj$fixed must be a list of numeric matrices (NAs ok).\n") 

  ## Set m and n
  if(is.null(modelObj$fixed$Z) || length(dim(modelObj$fixed$Z)) != 2) {
    msg = c(msg, "Something is wrong with fixed$Z. It should be a matrix.\n")
  }
  else {
    n = dim(modelObj$fixed$Z)[1]
    m = dim(modelObj$fixed$Z)[2]
  }

###########################
# Check that fixed and free are complete and consistent
###########################
  if(is.null(msg)) {  #if there is a problem with modes or presence then don't run further tests

  ## Set up element names that appear as fixed/free pair
  en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")
  
  ## Correct dimensions for reporting
  correct.dim1 = c(n,n,n,m,m,m,m,m)
  correct.dim2 = c(m,1,n,m,1,m,1,m)
  names(correct.dim1) = names(correct.dim2) = en

  fix = est = dim.fix = dim.est = nomatch = NULL

  for (elem in en) {

  ## Check for problems in the fixed/free pairs. Problems show up as TRUE 
    fix.null.flag <- is.null(modelObj$fixed[[elem]]) 
    est.null.flag <- is.null(modelObj$free[[elem]])
    dim.fix.flag = dim.est.flag = nomatch.flag = FALSE
    nonnum.flag = nonchar.flag = FALSE
 
    if (!fix.null.flag && !est.null.flag) {
    # check dim
      dim.fix.flag <- MARSScheckdims(elem, modelObj$fixed, n, m)
      dim.est.flag <- MARSScheckdims(elem, modelObj$free, n, m) 

    # check NAs match up; numeric
      if (!dim.fix.flag && !dim.est.flag) {
        fix.na <- which(is.na(modelObj$fixed[[elem]]))
        free.na <- which(is.na(modelObj$free[[elem]]))
        fix.num <- which(!is.na(modelObj$fixed[[elem]]))
        free.char <- which(!is.na(modelObj$free[[elem]]))
	    
        # Every NA in fixed must match up with a non-NA value in free
        # so no NAs can overlap and no non-NAs can overlap
        if(length(intersect(fix.na,free.na))!=0) nomatch.flag = TRUE
        if(length(intersect(fix.num,free.char))!=0) nomatch.flag = TRUE
      }
    }
 
    fix <- c(fix, fix.null.flag)
    est <- c(est, est.null.flag)
    dim.fix <- c(dim.fix, dim.fix.flag)
    dim.est <- c(dim.est, dim.est.flag) 
    nomatch <- c(nomatch, nomatch.flag)  
  }
  
  problem <- any(c(fix, est, dim.fix, dim.est, nomatch))

  if (problem) {
    if(any(fix)) {
      msg = c(msg, paste("Missing fixed", en[fix],"\n"))
    }
    if(any(est)) {
      msg = c(msg, paste("Missing free", en[est],"\n"))
    }
    if(any(dim.fix)) {
      msg = c(msg, paste("Dimension problem in fixed", en[dim.fix], "Dims should be", correct.dim1[dim.fix], "x", correct.dim2[dim.fix],"\n"))
    }
    if(any(dim.est)) {
      msg = c(msg, paste("Dimension problem in free", en[dim.est], "Dims should be", correct.dim1[dim.est], "x", correct.dim2[dim.est],"\n"))
    }
    if(any(nomatch)) {
      msg = c(msg, paste("Problem with the fixed/free pair for ", en[nomatch], ". The NAs and non-NAs don't match up correctly.\n", sep=""))
    }
  } 

###########################
# Check that none of the var-cov matrices have negative values on the diagonal
###########################
  en = c("R", "Q", "V0")
  neg = NULL
  for (elem in en) {
    neg.flag = FALSE
    if(any(takediag(modelObj$fixed[[elem]])<0,na.rm=TRUE)) neg.flag=TRUE
    neg = c(neg, neg.flag)
    }
  if(any(neg)) {
      msg = c(msg, paste("Negative values are on the diagonal of ", en[neg], ". Neg values are illegal on the diag of a var-cov matrix.\n", sep=""))
      }
      
###########################
# Check that the Q and R matrices have no zeros
###########################
  en = c("R", "Q")
  #this one catches if the user set 1 0 on the diagonal and set the others to NA (estimated)
  zer = NULL
  for (elem in en) {
    zero.flag = FALSE
    if(any(takediag(modelObj$fixed[[elem]])==0,na.rm=TRUE)) zero.flag=TRUE
    zer = c(zer, zero.flag)
    }
  if(any(zer)) { 
      msg = c(msg, paste("Zeros are on the diagonal of ", en[zer], ". In a MARSS model, zeros are illegal on the diag of the R or Q matrices.\n", sep=""))
      }

###########################
# Check that fixed V0, Q and R matrices are positive and definite
###########################
  en = c("R", "Q", "V0")
  zer = NULL
  for (elem in en) {
    zero.flag = FALSE
    #this if will not be triggered if all zero, but all zero Q or R will be caught above.
    if(is.fixed(modelObj$fixed[[elem]]) && !identical(unname(modelObj$fixed[[elem]]),array(0,dim=dim(modelObj$fixed[[elem]])))){
      tmp = try( eigen(modelObj$fixed[[elem]]), silent=TRUE )
      if(class(tmp)=="try-error") zero.flag=TRUE
      else if(any(tmp$values<=0)) zero.flag=TRUE
    }
    zer = c(zer, zero.flag)
  }
  if(any(zer)) msg = c(msg, paste("The fixed matrix ", en[zer], " is not positive-definite (and var-cov matrices must be).\n", sep=""))

        
###########################
# Check data and missing values consistency if data present
###########################
  if(!is.null(modelObj$data)) {

    if(!is.numeric(modelObj$data)) msg = paste(msg, "Data should be numeric. \n")
    if( (dim(modelObj$data)[1] != n) || length(dim(modelObj$data)) < 2)
      msg = c(msg, "Data is not a matrix or its dimensions do not match Z.\n")

    ## Check missing values matrix consistent with data 
    if(is.null(modelObj$miss.value) || is.null(modelObj$M)) {
      msg = c(msg, "Component miss.value or component M is missing.\n")
    }
    else {
      TT = dim(modelObj$data)[2]  
      for(i in 1:TT) {
        if( !isTRUE( all.equal(takediag(modelObj$M[,,i]), as.numeric(modelObj$data[,i] != modelObj$miss.value)) ) )
	  msg = c(msg, "Missing data matrix is inconsistent with data. Call marssm() to fix.\n")
      }
    }
  } # end if(!is.null(modelObj$data))

  } # end if(is.null(msg))

if(length(msg) == 0) return(TRUE)
else {
  msg=c("\nErrors were caught in is.marssm()\n", msg)
  return(msg)
}
}
