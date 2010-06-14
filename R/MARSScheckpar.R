#######################################################################################################
#   MARSScheckpar function
#   Utility function to check marssMLE params
#   Called by is.marssMLE(), MARSSsimulate()
#######################################################################################################
MARSScheckpar <- function(parList, n, m)
{    ## Set up element names
    en = c("Z", "A", "R", "B", "U", "Q", "x0", "V0")
    ## Correct dimensions for reporting
    correct.dim1 = c(n,n,n,m,m,m,m,m)
    correct.dim2 = c(m,1,n,m,1,m,1,m)
    names(correct.dim1) = names(correct.dim2) = en  
    param.null = dim.param = NULL
    msg = c()

    for (elem in en) {

      param.null.flag <- ( is.null(parList[[elem]]) || !is.finite(parList[[elem]]) )
      dim.param.flag = FALSE
 
      if (!param.null.flag) {
      	dim.param.flag <- MARSScheckdims(elem, parList, n, m) 
      }
 
      param.null <- c(param.null, param.null.flag)
      dim.param <- c(dim.param, dim.param.flag)  
    }
  
    problem <- any(c(param.null, dim.param))

    if (problem) { 
      if(any(param.null)) {
      	msg = c(msg, "MARSScheckpar: Missing or non-numeric parameter", en[param.null])
      }
      if(any(dim.param)) {
      	msg = c(msg, "MARSScheckpar: Dimension problem: parameter", en[dim.param], "Dims should be", correct.dim1[dim.param], "x", correct.dim2[dim.param])
      }
    }

    #check that Q and R are proper var-cov matrices
    if (!(all.equal(parList$Q, t(parList$Q)) && all(eigen(parList$Q)$values >= 0))) 
            msg = c(msg, "MARSScheckpar: Q is not a valid variance matrix")
    if (!(all.equal(parList$R, t(parList$R)) && all(eigen(parList$R)$values >= 0))) 
            msg = c(msg, "MARSScheckpar: R is not a valid variance matrix")

    if(length(msg) == 0) return(TRUE)
    else return(msg)
} 
