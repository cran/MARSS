describe.marssm = function(modelObj){
#This returns the structure of a model using text strings
fixed = modelObj$fixed
free = modelObj$free
m=dim(fixed$Z)[2]; n=dim(fixed$Z)[1]

constr.type = list(Q="none",R="none", B="none", U="none", A="none", Z="none", "x0"="none", "V0"="none")

  ############ Set the constraint type for U, A, and Z
  for(elem in c("U","A","x0")) {
    dimm = dim(free[[elem]])[1]
    while(constr.type[[elem]]=="none") {
      if( elem=="A" & is.design(fixed$Z) ){
        okai=rep(NA,m)
        for(i in 1:m){
            ai = fixed$A[as.logical(fixed$Z[,i])]
            okai[i] = (sum(ai==0,na.rm=TRUE)==1 & sum(is.na(ai))==(length(ai)-1))
            }
        if(all(okai)) { constr.type[[elem]] = "scaling"; break } 
        }
      if(is.fixed(fixed[[elem]]))  { 
        if(all(fixed[[elem]]==0))  { constr.type[[elem]] = "zero"; break }
        constr.type[[elem]] = "fixed"; break }
      if(length(free[[elem]])==1)  { constr.type[[elem]] = "scalar"; break }
      tmp = unique(as.character(free[[elem]])); tmp = tmp[!is.na(tmp)]
      if(length(tmp)==dimm) { constr.type[[elem]] = "unconstrained"; break }
      if(!any(is.na(free[[elem]])) & length(tmp)==1) { constr.type[[elem]] = "equal"; break }
      constr.type[[elem]] = "see fixed/free"
   }
   }

  ############ Check free$Q  and free$R
  # This part determines the form of free/fixed
  for(elem in c("Q","R","B","Z","V0")) {
     dimm = dim(free[[elem]])[1]
     dimc = dim(free[[elem]])[2]
    while(constr.type[[elem]]=="none") {
      if(elem=="Z" & is.design(fixed[[elem]])) { 
          Z.names = c()
          tmp = colnames(fixed$Q); X.names=c()
          if(is.null(tmp)) tmp=as.character(1:m)
          for(i in 1:length(tmp))
            X.names=c(X.names, strsplit(tmp[i],":")[[1]][1])
          for(i in 1:n) Z.names = c(Z.names, X.names[as.logical(fixed$Z[i,])])
          constr.type[[elem]] = Z.names
          break }
      if(is.fixed(fixed[[elem]])) { #it's fixed
        if(all(fixed[[elem]]==0))  { constr.type[[elem]] = "zero"; break }
        if(is.identity(fixed[[elem]])) constr.type[[elem]] = "identity"
        else constr.type[[elem]] = "fixed"
        break
        }
      if(length(free[[elem]])==1)  { constr.type[[elem]] = "scalar"; break }
      tmp = unique(as.character(free[[elem]])); tmp = tmp[!is.na(tmp)]
      if(length(tmp)==(dimm*dimc)) { constr.type[[elem]] = "unconstrained"; break }
      tmp.free = free[[elem]] 
      #need to deal with case where user uses "0" as a name
      if(isTRUE(any(free[[elem]]==0))) tmp.free = as.numeric(tmp.free) + 1 
      tmp.freefix = tmp.free
      tmp.free[is.na(tmp.free)]=0  #now all free values are non-zero and fixed values are zero
      tmp.freefix[is.na(free[[elem]])] = fixed[[elem]][is.na(free[[elem]])] #replace the NA values with fixed values
      if(is.equaltri(free[[elem]])) { constr.type[[elem]] = "equalvarcov"; break }  #no NAs in free in this case
      if(is.diagonal(tmp.freefix)) {
         tmp = unique(as.character(takediag(free[[elem]]))); tmp.no.na = tmp[!is.na(tmp)]
         if(length(tmp.no.na)==dimm) { constr.type[[elem]] = "diagonal and unequal"; break }
         if(length(tmp)==1 && !is.na(tmp)) { constr.type[[elem]] = "diagonal and equal"; break }
         if(length(tmp)==0) { constr.type[[elem]] = "error in describe.marssm. this case shouldn't be reached"; break }
         if(length(tmp.no.na)==0) { constr.type[[elem]] = "error in describe.marssm. this case shouldn't be reached"; break }
         if(length(tmp)>1 && !any(is.na(tmp))) { constr.type[[elem]] = paste("diagonal with ",length(tmp)," groups",sep="") ; break }
         if(any(is.na(tmp))) { constr.type[[elem]] = paste("diagonal with fixed elements and ",length(tmp.no.na)," groups",sep="") ; break }
         }
      tmp.unique = is.blockunconst(tmp.free, uniqueblocks=TRUE) & is.blockunconst(tmp.freefix, uniqueblocks=TRUE)
      tmp.uniqueornot = is.blockunconst(tmp.free, uniqueblocks=FALSE) & is.blockunconst(tmp.freefix, uniqueblocks=FALSE)
      if(tmp.unique) { constr.type[[elem]] = "unique block diagonal unconstrained"; break }
      if( tmp.uniqueornot && !tmp.unique ) { constr.type[[elem]] = "shared block diagonal unconstrained"; break }      
      tmp.unique = is.blockequaltri(tmp.free, uniqueblocks=TRUE) & is.blockunconst(tmp.freefix, uniqueblocks=TRUE)
      tmp.uniqueornot = is.blockequaltri(tmp.free, uniqueblocks=FALSE) & is.blockunconst(tmp.freefix, uniqueblocks=FALSE)
      if(tmp.unique) { constr.type[[elem]] = "unique block diagonal equalvarcov"; break }
      if(tmp.uniqueornot && !tmp.unique) { constr.type[[elem]] = "shared block diagonal equalvarcov"; break }
      constr.type[[elem]]="see fixed/free" #not assigned to one of the above cases
      }
    }

return(constr.type)
}
