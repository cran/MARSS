###########################################################################################################################
#   diag helper functions
###########################################################################################################################
takediag = function(x)
############# Function to take the diagonal; deals with R trying to think too much with diag()
{
if(length(x)==1) return(x)
else return(diag(x))
}

makediag = function(x,nrow=NA)
############# Function to make a diagonal matrix; deals with R trying to think too much with diag()
{
if(length(x)==1) 
{
if(is.na(nrow)) nrow=1
return(diag(c(x),nrow))
}
if((is.matrix(x) | is.array(x)))
   if(!(dim(x)[1]==1 | dim(x)[2]==1))  stop("Error in call to makediag; x is not vector")
if(is.na(nrow)) nrow=length(x)
return(diag(c(x),nrow))
}

############ The following functions are for testing the shapes of matrices
# these functions use as.character(x) to deal with the "feature" that in R (.01 + .14) == .15 is FALSE!!
is.equaltri = function(x) {
#warning this returns TRUE if x is 1x1
x = as.matrix(unname(x))
if(dim(x)[1]==1 & dim(x)[2]==1) return(TRUE)
#equal and non zero on diagonal; equal (zero ok) but different from diagonal on off-diagonal
if(any(is.na(x))) return(FALSE)
if(length(unique(as.character(x)))!=2) return(FALSE) #must be only 2 numbers in the matrix
if(is.diagonal(x)) return(TRUE)  #diagonal is special case of equaltri
#not diagonal
diagx = takediag(x)
tmp=table(as.character(diagx))
namestmp = names(tmp)
if(length(tmp)==1){
 if(length(unique(as.character(x)))!=2) return(FALSE) #not block diagonal
 if(length(unique(as.character(x)))==2) return(TRUE) 
}
return(FALSE) #length(tmp)!=1; must be only 1 number on diagonal
}

is.diagonal = function(x, na.rm=FALSE) {
#na.rm=TRUE means that NAs on the DIAGONAL are ignored
#non zero on diagonal; zero on off-diagonals
x=as.matrix(unname(x))
if(na.rm==FALSE && any(is.na(x))) return(FALSE)
nr = dim(x)[1]; nc = dim(x)[2];
if(nr != nc) return(FALSE)
diagx = takediag(x)
#ok if there are 0s on diagonal
#if(isTRUE( any(diagx==0) ) ) return(FALSE)
dimx = dim(x)[1]
if(length(x)==1) return(TRUE)
if(isTRUE( all(x[lower.tri(x)]==0) ) && isTRUE( all(x[upper.tri(x)]==0) ) ) return(TRUE) #diagonal
return(FALSE)
}

zero.row.col = function(x) {
#non zero on diagonal; zero on off-diagonals
x=as.matrix(unname(x))
if(!is.numeric(x)) return(FALSE)
diagx=takediag(x)
if(all(is.na(diagx))) return(FALSE)
if(!any(diagx==0)) return(FALSE)
whichzero = which(diagx==0)
if(!all(x[whichzero,]==0)) return(FALSE)
if(!all(x[,whichzero]==0)) return(FALSE)
return(TRUE)
}

is.identity = function(x) {
if(!is.matrix(x)) stop("argument must be a matrix")
if(!is.diagonal(x)) return(FALSE)
if(!all(takediag(x) %in% c(1))) return(FALSE)
return(TRUE)
}

is.blockdiag = function(x) {
x=as.matrix(unname(x))
if(any(is.na(x))) return(FALSE)
nr = dim(x)[1]; nc = dim(x)[2];
if(nr != nc) return(FALSE)
diagx = takediag(x)
if(any(diagx==0)) return(FALSE)

#special cases. 1. diagonal
if(is.diagonal(x)) return(TRUE)

#special cases. 2. all non-zero
if(all(x!=0)) return(TRUE)
 
tmpx = x
for(i in 1:nr){
   block = which(tmpx[1,]!=0)
   if(any(tmpx[block,block]==0) ) return(FALSE)
   if(!any(tmpx[1,]==0)) break
   notblock = which(tmpx[1,]==0)
   tmpx=tmpx[notblock,notblock]
   dim(tmpx) = c(length(notblock),length(notblock))
   }

return(TRUE)
}

is.blockequaltri = function(x, uniqueblocks=FALSE) {
#this looks for a block diagonal matrix with same number on diag and same (but different) number on off-diagonal 
#warning this returns true for a 1x1 matrix, diagonal matrix, and equaltri matrix
x=as.matrix(unname(x))
if(is.equaltri(x)) return(TRUE) #equaltri is a special case of block equaltri
if(!is.blockdiag(x)) return(FALSE) #blockequaltri is a special case of block diag
nr = dim(x)[1]
tmpx = x
trivals = c() #holder for the values
for(i in 1:nr){
   block = which(tmpx[1,]!=0)
   if(!is.equaltri(tmpx[block,block]) ) return(FALSE)
   trivals = c(trivals,tmpx[block[1],block[1]])
   if(length(block)>1) trivals = c(trivals,tmpx[block[1],block[2]]) 
   if(!any(tmpx[1,]==0)) break
   notblock = which(tmpx[1,]==0)
   tmpx=tmpx[notblock,notblock]
   dim(tmpx) = c(length(notblock),length(notblock))
   }
if(uniqueblocks==TRUE && any(table(as.character(trivals),exclude=c(NA,NaN,0))>1) )  return(FALSE) #all must be unique
return(TRUE) #got through the check without returning FALSE, so OK
}

is.blockunconst = function(x, uniqueblocks=FALSE) {
#this looks for a block diagonal matrix with each block an unconstrained matrix
x=as.matrix(unname(x))
if(!is.blockdiag(x)) return(FALSE) #blockunconst is a special case of block diag
if(uniqueblocks==TRUE && any(table(as.character(x),exclude=c(NA,NaN,0))>1) )  return(FALSE) #all must be unique
nr = dim(x)[1]
tmpx = x
for(i in 1:nr){
   block = which(tmpx[1,]!=0)
   tmp = unique(tmpx[block,block])
   dimblock = length(block)
   if(!length(tmp)==(dimblock*dimblock) ) return(FALSE)
   if(!any(tmpx[1,]==0)) break
   notblock = which(tmpx[1,]==0)
   tmpx=tmpx[notblock,notblock]
   dim(tmpx) = c(length(notblock),length(notblock))
   }
return(TRUE) #got through the check without returning FALSE, so OK
}

is.design = function(x) {
if(!is.matrix(x)) return(FALSE)
x=as.matrix(unname(x))
if(any(is.na(x)) || any(is.nan(x))) return(FALSE)
if(!is.numeric(x)) return(FALSE)
if(!all(x %in% c(1,0))) return(FALSE)
if(dim(x)[1]<dim(x)[2]) return(FALSE) #if fewer rows than columns then not design
tmp = rowSums(x)
if(!isTRUE(all.equal(tmp,rep(1,length(tmp))))) return(FALSE)
return(TRUE)
}

is.fixed = function(x) {
if(!is.matrix(x)) return(FALSE)
if(any(is.na(x))) return(FALSE)
if(!is.numeric(x)) return(FALSE)
return(TRUE)
}

vec = function(x) {
if(!is.matrix(x)) stop("arg must be a matrix")
return(matrix(x,length(x),1))
}

unvec = function(x,dim=NULL){
if(!is.vector(x) & !is.matrix(x)) stop("arg must be a vector or nx1 matrix)")
if(is.matrix(x))
  if(dim(x)[2]!=1) stop("if arg is a matrix it must be nx1 matrix)")
if(is.null(dim)) dim=c(length(x),1)
if(!is.vector(dim) & length(dim)!=2) stop("dim must be a vector of length 2: c(nrows,ncols)")
if(!is.numeric(dim)) stop("dim must be numeric")
if(!all(is.wholenumber(dim))) stop("dim must be a vector of 2 integers")
if(dim[1]*dim[2]!=length(x)) stop("num elements in arg greater than dim[1]*dim[2]")
return(array(x,dim=dim))
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5) {
  if(!is.numeric(x)) return(FALSE)
  test = abs(x - round(x)) < tol
  if(any(is.na(test))) return(FALSE)
  return(test)
  }

as.design = function(fixed,free) {
  if(!all(is.na(fixed)) & (!is.numeric(fixed) | !is.matrix(fixed))) stop("fixed argument must be a numeric matrix")
  if(!all(is.na(free)) & ( !is.character(free) | !is.matrix(free)) ) stop("free argument must be a character or NA matrix")
  f=vec(fixed); f[is.na(f)]=0
  tmp=table(free, exclude=c(NA,NaN))
  est.levels = names(tmp)
  numGroups = length(est.levels)      
  D.mat = matrix(0,dim(free)[1]*dim(free)[2],numGroups)   # matrix to allow shared growth rates (called F in my write-up)
  for(i in est.levels) D.mat[which(free==i),which(est.levels==i)] = 1
  if(all(is.na(free))) D.mat=array(0,dim=c(dim(free)[1]*dim(free)[2],1))
return(list(f=f, D=D.mat))
}

Imat = function(x) return(diag(1,x))

rwishart=function (nu, V) 
{
#function adapted from bayesm package
#author Peter Rossi, Graduate School of Business, University of Chicago
    m = nrow(V)
    df = (nu + nu - m + 1) - (nu - m + 1):nu
    if (m > 1) {
        T = diag(sqrt(rchisq(c(rep(1, m)), df)))
        T[lower.tri(T)] = rnorm((m * (m + 1)/2 - m))
    }else {
        T = sqrt(rchisq(1, df))
    }
    U = chol(V)
    C = t(T) %*% U
    return(crossprod(C))
}