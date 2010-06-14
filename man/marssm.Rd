\name{marssm}
\alias{marssm}
\alias{is.marssm}
\alias{as.marssm}
\title{ Model Objects }
\description{
These are model objects and utility functions for model objects in the package \code{\link{MARSS-package}}.
  \code{marssm()} creates multivariate autoregressive state space model objects. 
  \code{is.marssm()} ensures model consistency.
  \code{as.marssm()} attempts to coerce its argument to a marssm object.
}
\usage{
marssm(data = NULL, fixed, free, miss.value = -99)
is.marssm(modelObj)
as.marssm(wrapperObj)
}

\arguments{
  \item{data}{ An optional matrix (not dataframe), sites (rows) x years (columns), of observed population abundances.  If the algorithm is to be applied to log-abundance, the log transformation should be done before the data is passed in. }
  \item{fixed}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0, specifying which elements of each parameter are fixed. See Details.}
  \item{free}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0, specifying which elements of each parameter are to be estimated. See Details. }
  \item{miss.value}{ How are missing values represented in the data? }
  \item{modelObj}{ An object of class marssm. }
  \item{wrapperObj}{ An object of a wrapper class \code{\link{popWrap}}. }
}
\details{
A \code{marssm} object is an R representation of a MARSS model along with the data.  
    Data for \code{marssm()} consists of multivariate time series data in which time is across columns and the n observed time series are in the n different rows.

The MARSS model is
\describe{
\item{x(t+1) = B x(t) + U + w(t), where w(t) ~ MVN(0,Q)}{}
\item{y(t) = Z x(t) + A + v(t), where v(t) ~ MVN(0,R)}{}
\item{x(1) ~ MVN(x0, V0)}{} 
}
In the marssm object, the MARSS model is specified by fixed/free pairs of matrices
  for each parameter: B, U, Q, Z, A, R, x0, V0.  
    The dimensions for \code{fixed} and \code{free} matrices must be:
  \describe{
    \item{Z}{ n x m (m<=n)}
    \item{B}{ m x m }
    \item{U}{ m x 1 }
    \item{Q}{ m x m }
    \item{A}{ n x 1 }
    \item{R}{ n x n }
    \item{x0}{ m x 1 }
    \item{V0}{ m x m }
  }
 The matrices in \code{fixed} and \code{free} work as pairs to specify the fixed and free elements for each parameter. 
 In \code{fixed}, fixed elements must be numbers; values that are not fixed (i.e. are to be estimated) should be denoted NA.
  Elements in \code{free} will be interpreted as names for the free elements (even if they are numbers). 
  Identical elements within a parameter matrix will be constrained to have the same value. Non-free (i.e. fixed) values should 
  be denoted with NA, not 0, since the code will interpret 0 as "0" and assume that the user wants those parameters to be 
  coded with the name "0" and to be estimated.    See the manual (\code{show.doc(MARSS, manual)}) for many examples of MARSS models and their specification.
}
\value{
  An object of class "marssm".
  \item{data}{ Data supplied by user. }
  \item{fixed}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{free}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{M}{ An array of dim n x n x t (an n x n missing values matrix for each time point).  Each matrix is diagonal with 0 at the i,i value if the i-th value of y is missing, and 1 otherwise.} 
  \item{miss.value}{ Specifies missing value representation. Default is -99}
}
\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
}

\seealso{ \code{\link{popWrap}} } 
\examples{

  n = m = 5
  fixed = free = vector("list")

  free$Q=array(seq(1,m*m),dim=c(m,m)); fixed$Q=array(NA,dim=c(m,m))
  free$R=array(NA,dim=c(n,n)); diag(free$R)=1
  fixed$R=array(0,dim=c(n,n)); diag(fixed$R)=NA
  free$Z=array(NA,dim=c(m,m))
  fixed$Z=array(0,dim=c(m,m)); diag(fixed$Z)=1
  free$U=array(seq(1,m),dim=c(m,1)); fixed$U = array(NA,dim=c(m,1))
  fixed$A = matrix(NA,n,1); free$A = matrix(1:n,n,1)  
  free$B=array(NA,dim=c(m,m));
  fixed$B=array(0,dim=c(m,m)); diag(fixed$B)=1
  free$x0=array(seq(1,m),dim=c(m,1)); fixed$x0 = array(NA,dim=c(m,1));
  free$V0=array(NA,dim=c(m,m)); fixed$V0 = array(0,dim=c(m,m))

  m1 <- marssm(fixed=fixed, free=free)
  is.marssm(m1)

  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  aa = MARSS:::allowed$kem
  wrapperObj = popWrap(dat, allowed=aa, method="kem")
  modelObj = as.marssm(wrapperObj)

}

