\name{MARSSinits}
\alias{MARSSinits}
\title{ Initial Values for MLE }
\description{
  Sets up generic starting values for parameters for maximum-likelihood estimation algorithms that use an iterative maximization routine needing starting values.   Examples of such algorithms are the EM algorithm in \code{\link{MARSSkem}} and Newton methods in \code{\link{MARSSoptim}}.   This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSinits(modelObj, inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=5), method)
}
\arguments{
  \item{modelObj}{ An object of class \code{\link{marssm}}. \code{MARSSinits} uses three elements of the model object.
        \itemize{
    \item \code{data}  The data element is used to determine n, the dimension of the y in the MARSS model.  
    \item \code{fixed}  The fixed matrices are used to determine which (if any) model parameters are fixed. 
    \item \code{miss.value}  Used if a linear regression through the data is used to construct inits for x0.  In this case, the missing values are replaced by NA.  
    }
}
  \item{inits}{ A list of up to 8 values to construct starting matrices for each parameter in a MARSSmodel.}
  \item{method}{ The fitting algorithm to be used (e.g. kem) }
}
\details{
  Creates an \code{inits} parameter list for use by iterative maximization algorithms.  
  
  Defaults values for \code{inits} is supplied in \code{MARSSsettings.R}. The user can alter these and supply any of the following (m is the dim of X and n is the dim of Y in the MARSS model): 
\itemize{
    \item{elem=\code{A,U}} A numeric vector or matrix which will be constructed into \code{inits$elem} by the command \code{array(inits$elem),dim=c(n or m,1))}. If elem is fixed in the model, any \code{inits$elem} values will be overridden and replaced with the fixed value. Default is \code{array(0,dim=c(n or m,1))}. 
    \item{elem=\code{Q,R,B}} A numeric vector or matrix.  If length equals the length modelObj$fixed$elem then inits$elem will be constructed by \code{array(inits$elem),dim=dim(modelObj$fixed$elem))}. If length is 1 or equals m or n then inits$elem will be constructed into a diagonal matrix by the command \code{diag(inits$elem,m or n)}. If elem is fixed in the model, any \code{inits$elem} values will be overridden and replaced with the fixed value. Default is \code{diag(0.05,m or n)} for Q and R.  Default is \code{diag(1,m)} for B.
    \item{\code{x0}} If \code{inits$x0=-99}, then starting values for x0 are estimated by a linear regression 
    through the count data assuming A=0.  This will be a poor start if \code{inits$A} is not 0. If \code{inits$x0} is a numeric vector or matrix, \code{inits$x0} will be constructed by the command \code{array(inits$x0),dim=c(m,1))}. If x0 is fixed in the model, any \code{inits$x0} values will be overridden and replaced with the fixed value. Default is \code{inits$x0=-99}. 
    \item{\code{Z}} If Z is fixed in the model, \code{inits$Z} set to the fixed value. If Z is not fixed, then the user must supply \code{inits$Z}.  There is no default.
    \item{elem=\code{V0}} V0 is never estimated, so this is never used.
}
}
\value{
  A list with 8 matrices A, R, B, U, Q, x0, V0, Z, specifying initial values for parameters in a MARSS model.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{marssm}} \code{\link{MARSSkem}} \code{\link{MARSSoptim}} }


