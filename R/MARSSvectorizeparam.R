#######################################################################################################
#   MARSSvectorizeparam  function
#   Returns a vector of the ESTIMATED parameters or if vector passed in, that is put into list form for MLEobj$model
#######################################################################################################
MARSSvectorizeparam = function(MLEobj, parvec=NA) {
#This helper function
  #if parvec=NA) returns a vector version of all the estimated parameters (for use in say optim) from a mssm  model
  #if parvec is passed in) returns a marssMLE object with par fixed by parvec
  
  ## Order of elements is alphabetical
  en = c("A", "B", "Q", "R", "U", "V0", "x0", "Z")
  free = MLEobj$model$free
  fixed = MLEobj$model$fixed
  param = MLEobj$par

  ## If parvec = NA
  if(length(parvec)==1 && is.na(parvec[1])) {
    paramvector = NULL

    for(elem in en) {
      this.free = free[[elem]]
      this.param = param[[elem]]
      mat.names = unique(as.vector(this.free[!is.na(this.free)]))
      ## match() finds first appearance of each 
      matchvec = match(mat.names, this.free)
      tmp = this.param[matchvec]
      mat.names = paste(rep(elem, length(mat.names)), rep(".", length(mat.names)), mat.names, sep="")
      names(tmp) = mat.names
      paramvector = c(paramvector, tmp)
    }
    return(paramvector)
  } # end if parvec==NA


  ## If parvec passed in, use to set MLEobj$par	
  else { 
    parlen = 0; maxvec = NULL; par = NULL

    ## Check length(parvec) matches number of free params
    for(elem in en) {
      this.free = free[[elem]] 
      mx = length(unique(as.vector(this.free[!is.na(this.free)])))
      parlen = parlen + mx
      maxvec = c(maxvec, mx)
    }
    if(length(parvec) != parlen) stop("Length of param vector does not match # of free params")
    names(maxvec) = en
    
    ## Fill in values, matrix by matrix
    for(elem in en) {

      this.free = free[[elem]]
      this.fixed = fixed[[elem]]
      if(any(!is.na(this.free))) { ## if any free params
        mx = maxvec[[elem]]

        ## Params for this element
        elemvec = parvec[1:mx]
	if (is.null(names(elemvec))) {
	  num.names = unique(as.vector(this.free[!is.na(this.free)]))
        }
	else {
	  ## remove prefix
	  prefix = paste(elem, ".", sep="")
	  num.names = sub(prefix, "", names(elemvec))
	}
	## check name match
        if(!all(num.names %in% this.free)) 
	  stop(paste("parvec names don't match model$free names in parameter", elem))      
        ## Remove "used" values from parvec
        parvec = parvec[(mx+1):length(parvec)]

        ## Match names in matrix to param vector names  
        matchvec = match(this.free, num.names)   
        ## elemvec[matchvec] gets param values as vector
	#tmp = elemvec[matchvec]
	tmp = matrix(elemvec[matchvec], nrow=nrow(this.free))

	## Incorporate the fixed params
	fixedidx = which(!is.na(this.fixed))
	tmp[fixedidx] = this.fixed[fixedidx]
        tmp = list(tmp); names(tmp) = elem
        par = c(par, tmp)
      } #end if(any(!is.na(this.free)))

  else { ## use fixed matrix
	   par = c(par, fixed[elem])
      }

    } #end elem loop
    
    MLEobj$par = par
    return(MLEobj)
  } #end if parvec arg is passed in
}
