\name{MARSS}
\alias{MARSS}
\alias{MARSSsettings}
\title{ Interface MARSS Model Specification and Estimation }
\description{
  A top-level \code{\link{MARSS-package}} function to perform model specification and estimation for multivariate autoregressive state-space (MARSS) models.  To open the user guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}.  To open an overview page with package information, type \code{RShowDoc("index",package="MARSS")}.    
  
MARSS models take the form:
\describe{
\item{x(t+1) = B x(t) + U + w(t), where w(t) ~ MVN(0,Q)}{}
\item{y(t) = Z x(t) + A + v(t), where v(t) ~ MVN(0,R)}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}{} 
}
\code{MARSS} provides an interface to the base \code{\link{MARSS-package}} functions so that users do not need to directly construct \code{\link{marssm}} and \code{\link{marssMLE}} objects.
}
\usage{
MARSS(y, 
    inits=NULL,
    model=NULL,
    miss.value=NA,
    method = "kem",
    fit=TRUE, 
    silent = FALSE,
    control = NULL) 
}

\arguments{
The default settings for the optional arguments are set in \code{MARSSsettings.R} and are given below in the details section.  

  \item{y}{ A n x T matrix of n time series over T time steps. }
  \item{inits}{ List with up to 7 matrices A, R, B, U, Q, x0, V0, specifying initial values for parameters.  
    \itemize{
    \item \code{B }  Initial value(s) for B matrix (length 1 or m x m) 
    \item \code{U }  Initial value(s) for U matrix (length 1 or m x 1) 
    \item \code{Q }  Initial value(s) for process error variance(s) (length 1 or m x m) 
    \item \code{A }  Initial value(s) for observation bias (length 1 or n x 1) 
    \item \code{R }  Initial value(s) for non-process (observation) error variance(s) (length 1 or n x n) 
    \item \code{x0} Initial value(s) for hidden state(s) at time=1 or 0 (length 1 or m x 1) 
    \item \code{V0} Initial variance(s) for hidden state(s) at time=1 or 0 (length 1 or m x m) 
    }
  }
  \item{model}{ Model specification using parameter model descriptions. See Details. }
  \item{miss.value}{ How are missing values represented in the data? }
  \item{method}{ Estimation method. MARSS 1.0 allows method="kem" (see \code{\link{MARSSkem}} for more info) and "BFGS" (see \code{\link{MARSSoptim}} for more info). }
  \item{fit}{ TRUE/FALSE Whether to fit the model to the data. If FALSE, a marssMLE object with only the model is returned. }
  \item{silent}{ TRUE/FALSE Suppresses printing of full error messages, warnings, progress bars and convergence information. }
  \item{control}{ Estimation options for the maximization algorithm. The typically used control options for method="kem" are below but see  \code{\link{marssMLE}} for the full list of control options.  Note many of these are not allowed if method="BFGS"; see  \code{\link{MARSSoptim}} for the allowed control options for this method.
    \itemize{
      \item{\code{minit}}{   The minimum number of iterations to do in the maximization routine (if needed by method). If \code{method="kem"}, this is an easy way to up the iterations and see how your estimates are converging. (positive integer)} 
      \item{\code{maxit}}{   Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). }
      \item{\code{kf.x0}}{ Whether to set the prior at t=0 ("x00") or at t=1 ("x10"). The default is "x00". } 
      \item{\code{min.iter.conv.test}}{ Minimum iterations to run before testing convergence via the slope of the log parameter versus log iterations.} 
      \item{\code{conv.test.deltaT=9}}{ Number of iterations to use for the testing convergence via the slope of the log parameter versus log iterations. }
      \item{\code{conv.test.slope.tol}}{ The slope of the log parameter versus log iteration to use as the cut-off for convergence.  The default is 0.5 which is a bit high.  For final analyses, this should be set lower. }
      \item{\code{abstol}}{ Optional convergence tolerance for the maximization routine. If set to something other than NULL, then the absolute tolerance will be used as the convergence test.  Otherwise, the slope of the log of the parameters versus the log iteration is used as the test.} 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{   A positive integer.  If not 0, a record will be created during maximization iterations (what's recorded depends on the method of maximization). See \code{\link{MARSSkem}} and \code{\link{MARSSoptim}} for trace options and output details.}
      \item{\code{MCInit}}{   If TRUE, do a Monte Carlo search of the initial condition space. (T/F) } 
      \item{\code{numInits}}{   Number of random initial value draws if MCInit=TRUE (ignored otherwise). (positive integer)} 
      \item{\code{numInitSteps}}{   Number of EM iterations for each random initial value draw if MCInit=TRUE (ignored otherwise). (positive integer)}
      \item{\code{boundsInits}}{   Bounds on the uniform distributions from which initial values will be drawn if MCInit=TRUE (ignored otherwise). The Q (and analogously R) random matrices are created by a random draw from a wishart distribution where df=bounds[1] and S=diag(bounds[2],m). }
      \item{\code{silent}}{   1 or TRUE, Suppresses all printing including progress bars, error messages and convergence information. 0, Turns on all printing of progress bars, fitting information and error messages.  2, Prints a brief success/failure message.    }
  }
}
}
\details{
 MARSS provides an interface to the base \code{\link{MARSS-package}} functions and allows specification and fitting of MARSS models. The available estimation methods are maximum-likelihood via an EM algorithm (\code{method="kem"}) or via a quasi-Newton algorithm provided by function \code{optim} (\code{method="BFGS"}).  The function \code{MARSS()} allows the user to specify models using text strings for common classes of parameter matrices via the argument \code{model}.  It allows the user to specify fixed values for matrices by passing in numeric matrices in the \code{model} list.  Models with a combination of fixed values and shared estimated values can be specified using list matrices. See \code{\link{marssm}} or the user guide (reference and link below) for documentation and instructions on specifying list matrices for the model parameters.   

Valid model structures for \code{method="kem"} are below.  See the user guide (reference and link below) for details and type \code{allowed$kem} to see the allowed list specified in \code{MARSSsettings.R}.
\itemize{
    \item{\code{Z}} {A text string, "identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or "onestate", or a length n vector of factors specifying which of the m hidden state time series correspond to which of the n observation time series.   May be specified as a n x m list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric n x m matrix to use a custom fixed Z. "onestate" gives a n x 1 matrix of 1s. "identity","unconstrained", "diagonal and unequal", "diagonal and equal", and "equalvarcov" all specify m x m matrices.}
    \item{\code{B}} {"identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use custom fixed B, but in this case all the eigenvalues of B must fall in the unit circle.}
    \item{\code{U}, \code{A} and \code{x0}} {"unconstrained", "equal", "unequal" or "zero". May be specified as a m x 1 (or n x 1 for A) list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric m x 1 (or n x 1) matrix to use a custom fixed U (or A).  NAs can be put in this matrix to allow some elements to be fixed and others (the NAs) to be estimated.}
    \item{\code{Q}, \code{R} and \code{V0}} {"identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use a custom fixed matrix.  If the matrix is diagonal (off-diagonals all zeros), then NAs may appear on the diagonal to allow some diagonal elements to be fixed while other elements (the NAs) are estimated.}
    \item{special \code{A} option and comments} {The default setting for A is "scaling". This is used to treat A as an intercept where one A for each X (hidden state) is fixed at 0 and any other Ys associated with that X have an estimated A value.   Care must be taken when specifying A so that the model is not under-constrained and unsolveable model; i.e. at least one A element per X state needs to be fixed.}
  }       
  Valid model structures for \code{method="BFGS"} are the same as for \code{method="kem"}.  See  \code{\link{MARSSoptim}} for the allowed options for this method.

 The default estimation method, \code{method="kem"}, is the EM algorithm described in the user guide.  The default settings for the optional arguments are set via \code{MARSS:::alldefaults$kem} in \code{MARSSsettings.R}.  For this method, they are:
\itemize{
    \item{inits = list(B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=1)}
    \item{model = list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained")}
    \item{miss.value = NA}
    \item{control=list(minit=15, maxit=500, abstol=NULL, trace=0,
safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, MCInit=FALSE, numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,
boundsInits=list(B=c(0,1), U=c(-1,1), Q = c(sqrt(0.1),0.1,0.1),
Z=c(0,1), A=c(-1,1), R = c(sqrt(0.1),0.1,0.1) ) )}
  }
For \code{method="BFGS"}, type  \code{MARSS:::alldefaults$BFGS} to see the defaults.
      
 The likelihood surface for MARSS models can be  multimodal.  It is  recommended that for final analyses the ML estimates are checked by using the Monte Carlo initial conditions search using \code{MCInit=TRUE} in the \code{control} list.  This requires more computation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.
}

\value{
  An object of class marssMLE with the following components:
  \item{model}{ MARSS model specification (an object of class \code{\link{marssm}}). }
  \item{start}{ List with 8 matrices A, R, B, U, Q, Z, x0, V0, specifying initial values for parameters. } 
  \item{control}{ A list of estimation options, as specified by arguments \code{control}. }
  \item{method}{ Estimation method. }
  If \code{fit=TRUE}, the following are also added to the marssMLE object.
  If \code{fit=FALSE}, an marssMLE object ready for fitting via the specified \code{method} is returned.
  \item{par}{ A list with 8 matrices of estimated (+ fixed) parameter values Z, A, R, B, U, Q, x0, V0. }
  \item{kf}{ A list containing Kalman filter/smoother output from \code{\link{MARSSkf}}.}
  \item{numIter}{ Number of iterations required for convergence. }
  \item{convergence}{ Convergence status. 0 means converged successfully.  Anything else is a warning or error. 2 means the MLEobj has an error; the MLEobj is returned so you can debug it.  The other numbers are errors during fitting.  The error code depends on the fitting method.  See \code{\link{MARSSkem}} and \code{\link{MARSSoptim}}. }
  \item{logLik}{ Log-likelihood. }
  \item{AIC}{ Akaike's Information Criterion. }
  \item{AICc}{ Sample size corrected AIC. }
}
\references{ 
The user guide:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.          

Holmes, E. E. (2010).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Link on the CRAN page under EMDerivation.pdf
}
\author{ 
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, kellie(dot)wills(at)noaa(dot)gov
}
\seealso{
\code{\link{marssm}} \code{\link{marssMLE}} \code{\link{MARSSkem}} \code{\link{MARSSoptim}} \code{\link{MARSS-package}}
}
\examples{
#harborSealWA is a n=5 matrix of logged population counts
dat = t(harborSealWA)
dat = dat[2:nrow(dat),] #remove the year row
#fit a model with 1 hidden state and 5 observation time series
kemfit = MARSS(dat, model=list(Z=factor(c(1,1,1,1,1)), 
   R="diagonal and equal"))
kemfit$model #This gives a description of the model
print(kemfit$model) # same as kemfit$model
summary(kemfit$model) #This shows the model structure

#fit the model using quasi-Newton algorithm
\dontrun{ #takes a long time
bfgsfit = MARSS(dat, model=list(Z=factor(c(1,1,1,1,1)), 
   R="diagonal and unequal"), method="BFGS")
}

#add CIs to a marssMLE object
#default uses an estimated Hessian matrix
kem.with.hess.CIs = MARSSparamCIs(kemfit)   
kem.with.hess.CIs #print with se's and CIs

#estimate CIs using a parametric bootstrap
tmp=kemfit
tmp$control$abstol=0.01
kem.with.boot.CIs = MARSSparamCIs(tmp, 
   method="parametric", nboot=10)   
#nboot should be more like 1000, but set low for example sake
kem.with.boot.CIs #print with se's, CIs, and bias est

#fit a model with 5 hidden states (default)
kemfit = MARSS(dat, silent=TRUE)  #suppress printing
kemfit #print information on the marssMLE object

#fit a model with 5 correlated hidden states 
# with one variance and one  covariance using the delta logLik convergence test
# by setting abstol.  This is fast, but does not actually test convergence (i.e. that 
# log param value versus log iteration is flat
kemfit = MARSS(dat, model=list(Q="equalvarcov"), control=list(abstol=0.1))

#fit a model with 5 correlated hidden states
#abstol set rather huge and many of the Q's are not converged based on log-log test  
kemfit = MARSS(dat, model=list(Q="unconstrained"), control=list(abstol=0.1))

#fit a model with 5 independent hidden states 
#where each observation time series is independent
#the hidden trajectories 1-3 & 4-5 share their U parameter
kemfit = MARSS(dat, model=list(U=matrix(c("N","N","N","S","S"),5,1)))

#same model, but with fixed independent observation errors
#and the 4th and 5th x processes are forced to have a U=0
#Notice how a list matrix is used to combine fixed and estimated elements
#all parameters can be specified in this way using list matrices
kemfit = MARSS(dat, model=list(U=matrix(list("N","N","N",0,0),5,1), 
   R=diag(0.01,5)),control=list(minit=100))

#fit a model with 2 hidden states (north and south) 
#where observation time series 1-3 are north and 4-5 are south
#Make the hidden state process independent with same process var
#Make the observation errors different but independent
#Make the growth parameters (U) the same
kemfit = MARSS(dat, model=list(Z=factor(c("N","N","N","S","S")),
   Q="diagonal and equal",R="diagonal and unequal",U="equal"),
   control=list(minit=100))

#print the model followed by the marssMLE object
kemfit$model
kemfit

#simulate some new data from our fitted model
sim.data=MARSSsimulate(kemfit$par, nsim=10, tSteps=100)

\dontrun{
#Compute bootstrap AIC for the model; this takes a long, long time
kemfit.with.AICb = MARSSaic(kemfit, output = "AICbp")
kemfit.with.AICb
}
}

