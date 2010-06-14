\name{MARSS}
\alias{MARSS}
\alias{MARSSsettings}
\title{ Interface MARSS Model Specification and Estimation }
\description{
  A top-level \code{\link{MARSS-package}} function to perform model specification and estimation for multivariate autoregressive state-space (MARSS) models:
\describe{
\item{x(t+1) = B x(t) + U + w(t), where w(t) ~ MVN(0,Q)}{}
\item{y(t) = Z x(t) + A + v(t), where v(t) ~ MVN(0,R)}{}
\item{x(1) ~ MVN(x0, V0)}{} 
}
\code{MARSS} provides an interface to the base \code{\link{MARSS-package}} functions so that users do not need to directly construct \code{\link{marssm}} and \code{\link{marssMLE}} objects.
}
\usage{
MARSS(y, 
    inits=NULL,
    constraint=NULL,
    fixed=NULL, free=NULL, 
    miss.value=-99,
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
    \item \code{x0} Initial value(s) for hidden state(s) at time=1 (length 1 or m x 1) 
    \item \code{V0} Initial variance(s) for hidden state(s) at time=1 (length 1 or m x m) 
    }
  }
  \item{constraint}{ Model specification using parameter constraint descriptions. See Details. }
  \item{fixed}{ Optional model specification using matrices of fixed and free parameters. See Details. }
  \item{free}{ Optional model specification using matrices of fixed and free parameters. See Details. }
  \item{miss.value}{ How are missing values represented in the data? }
  \item{method}{ Estimation method. MARSS 1.0 allows method="kem" and "BFGS". }
  \item{fit}{ TRUE/FALSE Whether to fit the model to the data. If FALSE, a marssMLE object with only the model is returned. }
  \item{silent}{ TRUE/FALSE Suppresses printing of full error messages, warnings, progress bars and convergence information. }
  \item{control}{ Estimation options for the maximization algorithm. The control options for method="kem" are
    \itemize{
      \item{\code{minit}}{   The minimum number of iterations to do in the maximization routine (if needed by method). If \code{method="kem"}, this is an easy way to up the iterations and see how your estimates are converging. (positive integer)} 
      \item{\code{maxit}}{   Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). } 
      \item{\code{abstol}}{   Convergence tolerance for the maximization routine. (default is 0.01 which is a bit high.) } 
      \item{\code{iter.V0}}{   The value of V0 to be used in place of 0 when x0 is treated as fixed and V0=0. See manual for discussion of initial state variance. (default is 10 which works well)}
      \item{\code{trace}}{   A positive integer.  If not 0, a record will be created during maximization iterations (what's recorded depends on the method of maximization). }
      \item{\code{MCInit}}{   If TRUE, do a Monte Carlo search of the initial condition space. (T/F) } 
      \item{\code{numInits}}{   Number of random initial value draws if MCInit=TRUE (ignored otherwise). (positive integer)} 
      \item{\code{numInitSteps}}{   Number of EM iterations for each random initial value draw if MCInit=TRUE (ignored otherwise). (positive integer)}
      \item{\code{boundsInits}}{   Bounds on the uniform distributions from which initial values will be drawn if MCInit=TRUE (ignored otherwise). }
      \item{\code{silent}}{   1 or TRUE, Suppresses all printing including progress bars, error messages and convergence information. 0, Turns on all printing of progress bars, fitting information and error messages.  2, Prints a brief success/failure message.    }
  }
}
}
\details{
 MARSS provides an interface to the base \code{\link{MARSS-package}} functions and allows specification and fitting of MARSS models. In \code{\link{MARSS-package}} 1.0, the available estimation methods are maximum-likelihood via a Kalman-EM algorithm (\code{method="kem"}) or via a quasi-Newton algorithm provided by function \code{optim} (\code{method="BFGS"}).  The function \code{MARSS()} allows the user to specify models using text strings for common classes of parameter matrices via the argument \code{constraint}.  It allows the user to specify fixed values for matrices by passing in numeric matrices in the \code{constraint} list.  If the model classes available via the \code{constraint} strings are not sufficient, MARSS also allows specification using matrix pairs specified with argument \code{fixed} and \code{free}. If \code{fixed/free} matrices are specified for some parameters, these will override any constraints for those parameters.  See \code{\link{marssm}} or the manual (\code{show.doc(MARSS, manual)}) for documentation and instructions on specifying fixed and free matrices.   

Valid constraints for \code{method="kem"} are below.  See the manual (\code{show.doc(MARSS, manual)}) for details and type \code{allowed$kem} to see the allowed list specified in \code{MARSSsettings.R}.
\itemize{
    \item{\code{Z}} {"identity" or a vector of factors specifying which of the m hidden state time series correspond to which of the n observation time series.   May also be specified as a numeric n x m matrix to use a custom fixed Z.}
    \item{\code{B}} {"identity" or a vector of factors specifying shared diagonal elements.   May also be specified as a numeric m x m matrix to use custom fixed B, but in this case all the eigenvalues of B must fall in the unit circle.}
    \item{\code{U}} {"unconstrained", "equal", "unequal" or "zero". May also be a vector of factors specifying shared u terms.   May also be specified as a numeric m x 1 matrix to use a custom fixed U.  NAs can be put in this matrix to allow some elements to be fixed and others (the NAs) to be estimated.}
    \item{\code{Q}} {"unconstrained", "diagonal and unequal", "diagonal and equal", or "equalvarcov". May also be a vector of factors specifying shared diagonal elements.  May also be specified as a numeric m x m matrix to use a custom fixed Q.  If the matrix is diagonal (off-diagonals all zeros), then NAs may appear on the diagonal to allow some diagonal elements to be fixed while other elements (the NAs) are estimated.}
    \item{\code{A}} {"scaling" This treats A as an intercept or "zero" which sets A to a fixed value of all zeros.   May also be specified as a numeric n x 1 matrix to use a custom fixed A.  NAs can be put in this matrix to allow some elements to be fixed and others (the NAs) to be estimated.  Note all NAs in A would mean all A elements are estimated and would typically result in an under-constrained and unsolveable model; at least one A element per X state needs to be fixed.}
    \item{\code{R}} {"unconstrained", "diagonal and unequal", "diagonal and equal", or "equalvarcov". May also be a vector of factors specifying shared diagonal elements.  May also be specified as a numeric n x n matrix to use a custom fixed R.  If the matrix is diagonal (off-diagonals all zeros), then NAs may appear on the diagonal to allow some diagonal elements to be fixed while other elements (the NAs) are estimated.}
    \item{\code{x0}} {"unconstrained", "equal", "unequal" or "zero". May also be a vector of factors specifying shared initial state values.  May also be specified as a numeric m x 1 matrix to use a custom fixed x0. NAs can be put in this matrix to allow some elements to be fixed and others (the NAs) to be estimated.}
  }       
  Valid constraints for \code{method="BFGS"} are the same as for \code{method="kem"} except that the Q and R matrices must be diagonal if estimated.  Thus \code{"unconstrained"} and \code{"equalvarcov"} are not allowed for Q or R.  Type \code{allowed$BFGS} to see the allowed list specified in \code{\link{MARSSsettings}}.

 The default estimation method, \code{method="kem"}, is the "Kalman-EM" algorithm described in the manual.  The default settings for the optional arguments are set via \code{alldefaults$kem} in \code{MARSSsettings}.  For this method, they are:
\itemize{
    \item{inits = list(B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=10)}
    \item{constraint = list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained")}
    \item{miss.value = -99}
    \item{control=list(minit=1, maxit=5000, abstol=0.01, iter.V0=10, trace=0, safe=FALSE, MCInit=FALSE, numInits = 500, numInitSteps = 10, boundsInits=list(B=c(0,1), U=c(-1,1), 
      logQ = c(log(1.0e-05),log(1.0)), Z=c(0,1), A=c(-1,1), logR = c(log(1.0e-05),log(1.0))))}
  }
For \code{method="BFGS"}, type  \code{alldefaults$BFGS} to see the defaults.
      
 The likelihood surface for MARSS models can be highly multimodal.  It is  recommended that for final analyses the ML estimates are checked by using the Monte Carlo initial conditions search using \code{MCInit=TRUE} in the \code{control} list.  This requires more computation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.
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
The user manual:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   type \code{show.doc(MARSS, manual)} to see.          
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA.  

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
kemfit = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1)), 
   R="diagonal and unequal"))
kemfit$model #This gives a description of the model
print(kemfit$model) # same as kemfit$model
summary(kemfit$model) #This shows the model structure
show(kemfit$model)    #shows the raw object

#The fit above throws a warning about convergence
#Determine which variance has not yet converged
find.degenerate(kemfit)

\dontrun{
#Increase the number of iterations to ensure convergence
kemfit100 = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1)), 
   R="diagonal and unequal"), control=list(minit=100))
find.degenerate(kemfit100)  #now check convergence of variances

#fit the model using quasi-Newton algorithm
bfgsfit = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1)), 
   R="diagonal and unequal"), method="BFGS")

#add CIs to a marssMLE object
#default uses an estimated Hessian matrix
kem.with.hess.CIs = MARSSparamCIs(kemfit100)   
kem.with.hess.CIs #print with se's and CIs
#estimate CIs using a parametric bootstrap
kem.with.boot.CIs = MARSSparamCIs(kemfit100, 
   method="parametric", nboot=10)   
#nboot should be more like 1000, but set low for example sake
kem.with.boot.CIs #print with se's, CIs, and bias est

#fit a model with 5 hidden states (default)
kemfit = MARSS(dat, silent=TRUE)  #suppress printing
kemfit #print information on the marssMLE object
show(kemfit) #look at the raw marssMLE object

#fit a model with 5 correlated hidden states  
kemfit = MARSS(dat, constraint=list(Q="unconstrained"))

#fit a model with 5 correlated hidden states 
# with one variance and one  covariance
kemfit = MARSS(dat, constraint=list(Q="equalvarcov"))

#fit a model with 5 independent hidden states 
#where each observation time series is independent
#the hidden trajectories 1-3 & 4-5 share their U parameter
kemfit = MARSS(dat, constraint=list(U=factor(c("N","N","N","S","S"))))

#same model, but with fixed independent observation errors
kemfit = MARSS(dat, constraint=list(U=factor(c("N","N","N","S","S")), 
   R=diag(0.01,5)),control=list(minit=100))

#fit a model with 2 hidden states (north and south) 
#where observation time series 1-3 are north and 4-5 are south
#Make the hidden state process independent with same process var
#Make the observation errors different but independent
#Make the growth parameters (U) the same
kemfit = MARSS(dat, constraint=list(Z=factor(c("N","N","N","S","S")),
   Q="diagonal and equal",R="diagonal and unequal",U="equal"),
   control=list(minit=100))

#print the model followed by the marssMLE object
kemfit$model
kemfit

#simulate some new data from our fitted model
sim.data=MARSSsimulate(kemfit$par, nsim=10, tSteps=100)

#Compute bootstrap AIC for the model; this takes a long time
kemfit.with.AICb = MARSSaic(kemfit, output = "AICbp")
kemfit.with.AICb
}
}

