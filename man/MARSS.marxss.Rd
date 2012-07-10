\name{MARSS.marxss}
\alias{MARSS.marxss}
\title{ Multivariate AR-1 State-space Model with Inputs }
\description{
The MARXSS (MAR-1 with eXogenous variables) model is
\describe{
\item{x(t+1) = B(t) x(t) + u(t) + C(t)c(t) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + a(t) + D(t)d(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(0) ~ MVN(x0, V0)}{} 
}
}
\usage{
MARSS.marxss(MARSS.inputs)
}

\arguments{
  \item{MARSS.inputs}{ A list of the inputs to a MARSS() call. \code{MARSS.marxss()} uses the model argument.}
}
\details{
See the Quick Start guide (\code{RShowDoc("Quick_Start",package="MARSS")}) or the User Guide (\code{RShowDoc("UserGuide",package="MARSS")}) for examples.

\code{\link{MARSS}()} is called with the argument \code{model} which is a list.  The elements in the list specify the structure for the B, u, C, c, Q, Z, a, D, d, R, x0, and V0 in the MARXSS model (above).  The list elements can have the following values:
\itemize{
    \item{\code{Z}} {A text string, "identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or "onestate", or a length n vector of factors specifying which of the m hidden state time series correspond to which of the n observation time series.   May be specified as a n x m list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric n x m matrix to use a custom fixed Z. "onestate" gives a n x 1 matrix of 1s. "identity","unconstrained", "diagonal and unequal", "diagonal and equal", and "equalvarcov" all specify n x n matrices.}
    \item{\code{B}} {"identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use custom fixed B, but in this case all the eigenvalues of B must fall in the unit circle.}
    \item{\code{U}, \code{A} and \code{x0}} {"unconstrained", "equal", "unequal" or "zero". May be specified as a m x 1 (or n x 1 for A) list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric m x 1 (or n x 1) matrix to use a custom fixed U (or A).}
    \item{\code{Q}, \code{R} and \code{V0}} {"identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use a custom fixed matrix.}
    \item{special \code{A} option and comments} {The default setting for A is "scaling". This is used to treat A as an intercept where one A for each X (hidden state) is fixed at 0 and any other Ys associated with that X have an estimated A value.   Care must be taken when specifying A so that the model is not under-constrained and unsolveable model; i.e. at least one A element per X state needs to be fixed.}
    \item{\code{D} and \code{C}} {"identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use custom fixed values. Must have n rows (D) or m rows (C).}
   \item{\code{d} and \code{c}} {Numeric matrix of inputs.  No missing values allowed.  Must have 1 column or the same number of columns as the data, y.  The numbers of rows in d must be the same as number of columns in D; similarly for c and C.}
  }       
Z, A, D, R, B, U, C and Q may be time-varying.  If time-varying, then text shortcuts cannot be used.  Enter as an array with the 3rd dimension being time.   Time dimension must be 1 or equal to the number of time-steps in the data.    See Quick Start guide (\code{RShowDoc("Quick_Start",package="MARSS")}) or the User Guide (\code{RShowDoc("UserGuide",package="MARSS")}) for  examples.Valid model structures for \code{method="BFGS"} are the same as for \code{method="kem"}.  See  \code{\link{MARSSoptim}} for the allowed options for this method.

 The default estimation method, \code{method="kem"}, is the EM algorithm described in the user guide.  The default settings for the control and inits arguments are set via \code{MARSS:::alldefaults$kem} in \code{MARSSsettings.R}.  The defaults for the model argument are set in \code{MARSS.marxss.R} For this method, they are:
\itemize{
    \item{inits = list(Z=1, B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=0)}
    \item{model = list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero", C="zero",D="zero",c=matrix(0,0,1), d=matrix(0,0,1), tinitx=0, diffuse=FALSE)}
    \item{miss.value = NA}
    \item{control=list(minit=15, maxit=500, abstol=NULL, trace=0,
safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, MCInit=FALSE, numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,
boundsInits=list(B=c(0,1), U=c(-1,1), Q = c(sqrt(0.1),0.1,0.1),
Z=c(0,1), A=c(-1,1), R = c(sqrt(0.1),0.1,0.1) ) )}
  }
For \code{method="BFGS"}, type  \code{MARSS:::alldefaults$BFGS} to see the defaults.


}
\value{
  An object of class "marssm".
  \item{data}{ Data supplied by user. }
  \item{fixed}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{free}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{X.names}{ A m length vector of names for the X states.} 
  \item{miss.value}{ Specifies missing value representation. Default is NA}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  Kellie Wills, NOAA, Seattle, USA.  
  }

\seealso{ \code{\link{MARSS}} } 
