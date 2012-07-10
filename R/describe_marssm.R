describe.marssm = function(modelObj){
#This returns the structure of a model using text strings
fixed = modelObj$fixed
free = modelObj$free
m=dim(fixed$x0)[1]; n=dim(modelObj$data)[1]
model.elem=names(fixed)
constr.type = list(Q="none",R="none", B="none", U="none", A="none", Z="none", "x0"="none", "V0"="none")

for(elem in model.elem){
  TT.max = max( dim(fixed[[elem]])[3], dim(free[[elem]])[3])
  if(TT.max > 1) constr.type[[elem]]="time-varying"
}  
  # For everything below, elem is time constant however tests are written to take a 3D matrix and it is not req that t=1
    
  ############ Set the model structure for U, A, and x0
  for(elem in c("U","A","x0")) {
    dimm = dim(free[[elem]])[1]
    while(constr.type[[elem]]=="none") {
      if(is.fixed(free[[elem]]))  { 
        if(all(fixed[[elem]]==0))  { constr.type[[elem]] = "fixed and zero"; break }
        constr.type[[elem]] = "fixed"; break }
      if(length(free[[elem]])==1)  { constr.type[[elem]] = "scalar"; break }
      if(all(apply(free[[elem]],3,is.diagonal))) { constr.type[[elem]] = "unconstrained"; break }
      if(length(unique(as.vector(free[[elem]])))==1 & dim(free[[elem]])[2]==1) { constr.type[[elem]] = "equal"; break }
      constr.type[[elem]] = "use summary()"
   }
   }

  ############ Check the other matrices
  # This part determines the form of free/fixed
  for(elem in c("Q","R","B","Z","V0")) {
    if(elem == "Z"){
      dimm = dim(free$A)[1]
      dimc = dim(free$U)[1]
      }else{ dimm=dimc=sqrt(dim(free[[elem]])[1]) }
    while(identical(constr.type[[elem]],"none")) {
      #Z is not time varying and is a design matrix; note test above would have ensured Z is time constant
      if( elem=="Z" & is.fixed(free[[elem]]) & dim(fixed[[elem]])[3]==1 & all(apply(fixed[[elem]],3,is.design,dim=c(dimm,dimc)))) { 
          Z.names = c()
          tmp = colnames(fixed$Q); X.names=c()
          if(is.null(tmp)) tmp=as.character(1:m)
          for(i in 1:length(tmp))
            X.names=c(X.names, strsplit(tmp[i],":")[[1]][1])
          for(i in 1:n) Z.names = c(Z.names, X.names[as.logical(fixed$Z[i,,1])]) #fixed$Z will be time constant
          constr.type[[elem]] = Z.names
          break } #break out of the while for elem
      if(is.fixed(free[[elem]])) { #it's fixed
        if(all(fixed[[elem]]==0))  { constr.type[[elem]] = "fixed and zero"; break }
        if( all(apply(fixed[[elem]],3,is.identity,dim=c(dimm,dimc)))) constr.type[[elem]] = "identity"
        else constr.type[[elem]] = "fixed"
        break #break out of while
        }
      if(length(free[[elem]])==1)  { constr.type[[elem]] = "scalar"; break }
      if(all(apply(free[[elem]],3,is.identity))) { constr.type[[elem]] = "unconstrained"; break }
      tmp.free = free[[elem]] 

      #fixed.free.to.formula requires that 3D mats have dim3=1, which will be true here since check above that time constant
      #creates a list matrix version of model
      tmp.mat = fixed.free.to.formula( fixed[[elem]],free[[elem]],c(dimm,dimc) )  
      if(elem %in% c("Q","R","V0")){
        dimm=sqrt(dim(fixed[[elem]])[1])
        if(
          all(fixed[[elem]]==0) &
          all(unlist(tmp.mat[upper.tri(tmp.mat)])==unlist(t(tmp.mat)[upper.tri(tmp.mat)])) & 
          is.design(free[[elem]]) & 
          length(unique(unlist(tmp.mat[upper.tri(tmp.mat)])))==dimm*(dimm-1)/2 &
          length(unique(unlist(diag(tmp.mat))))==dimm )
             { constr.type[[elem]] = "unconstrained"; break } 
      }
      if(is.diagonal(tmp.mat)) {
         tmp.diag = unique(diag(tmp.mat))  #it's a list
         tmp.diag.est=tmp.diag[sapply(tmp.diag,is.character)]
         tmp.diag.fixed=tmp.diag[sapply(tmp.diag,is.numeric)]
         if(length(tmp.diag.est)==dimm) { constr.type[[elem]] = "diagonal and unequal"; break }
         if(length(tmp.diag.est)==1 & length(tmp.diag.fixed)==0) { constr.type[[elem]] = "diagonal and equal"; break }
         if(length(tmp.diag.est)>1 & length(tmp.diag.fixed)==0) { 
            constr.type[[elem]] = paste("diagonal with ",length(tmp.diag.est)," groups",sep="") ; break }
         if(length(tmp.diag.est)>1 & length(tmp.diag.fixed)>0) { 
            constr.type[[elem]] = paste("diagonal with fixed elements and ",length(tmp.diag.est)," groups",sep="") ; break }
         }
      if(is.equaltri(tmp.mat)) { 
         if(elem %in% c("Q","R","V0")){
           constr.type[[elem]] = "one variance value and covariance value"; break 
         }else{ constr.type[[elem]] = "one diagonal value and one off-diagonal value"; break }
      }  

      tmp.unique = is.blockunconst(tmp.mat, uniqueblocks=TRUE) & is.blockunconst(tmp.mat, uniqueblocks=TRUE)
      if(tmp.unique) { constr.type[[elem]] = "unique block diagonal unconstrained"; break }
      tmp.uniqueornot = is.blockunconst(tmp.free, uniqueblocks=FALSE) & is.blockunconst(tmp.mat, uniqueblocks=FALSE)
      if( tmp.uniqueornot && !tmp.unique ) { constr.type[[elem]] = "shared block diagonal unconstrained"; break }      
      tmp.unique = is.blockequaltri(tmp.mat, uniqueblocks=TRUE) & is.blockunconst(tmp.mat, uniqueblocks=TRUE)
      if(tmp.unique) { constr.type[[elem]] = "unique block diagonal equalvarcov"; break }
      tmp.uniqueornot = is.blockequaltri(tmp.mat, uniqueblocks=FALSE) & is.blockunconst(tmp.mat, uniqueblocks=FALSE)
      if(tmp.uniqueornot && !tmp.unique) { constr.type[[elem]] = "shared block diagonal equalvarcov"; break }
      constr.type[[elem]]="see summary()" #not assigned to one of the above cases
      }
    } #"Q","R","B","Z","V0"

return(constr.type)
}
