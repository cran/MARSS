\name{marssMLE}
\alias{marssMLE}
\alias{is.marssMLE}
\title{ Maximum Likelihood MARSS Estimation Object }
\description{
  An object in the  \code{\link{MARSS-package}} that has all the elements needed for maximum-likelihood estimation of multivariate autoregressive state-space model: the data, model, estimation methods, and any control options needed for the method.  If the model has been fit and parameters estimated, the object will also have the MLE parameters.  Other functions add other elements to the marssMLE object, such as CIs, s.e.'s, AICb, and hessian. 
}
\usage{
is.marssMLE(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. See Details.}
}
\details{
  The is.marssMLE() function checks components \code{model}, \code{start} and \code{control}, which must be present for estimation by functions e.g. \code{\link{MARSSkem}}. Components returned from estimation must include at least \code{method}, \code{par}, \code{kf}, \code{numIter}, \code{convergence} and \code{logLik}. Additional components (e.g. AIC) may be returned, as described in function help files.
\describe{
  \item{\code{model}}{ MARSS model specification (an object of class \code{\link{marssm}}). }
  \item{\code{start}}{ List with 7 matrices A, R, B, U, Q, x0, V0, specifying initial values for parameters to be used (if needed) by the maximization algorithm. 
    \describe{
    \item{\code{B}}{ Initial value(s) for B matrix (m x m). }
    \item{\code{U}}{ Initial value(s) for U matrix (m x 1). }
    \item{\code{Q}}{ Initial value(s) for Q variance-covariance matrix (m x m). }
    \item{\code{A}}{ Initial value(s) for A matrix (n x 1). }
    \item{\code{R}}{ Initial value(s) for R variance-covariance matrix (n x n). }
    \item{\code{x0}}{ Initial value(s) for initial state vector (m x 1). }
    \item{\code{V0}}{ Initial variance(s) for initial state variance (m x m). }
    }
  }
  \item{\code{control}}{ A list specifying estimation options.  The following options are needed by \code{\link{MARSSkem}}.  Other control options
  can be set if needed for other estimation methods, e.g. the control options listed for \code{\link{optim}} for use with \code{\link{MARSSoptim}}.  The default values for control options are set in \code{alldefaults[[method]]} which is specified in \code{MARSSsettings.R}.
    \describe{
      \item{\code{minit}}{ The minimum number of iterations to do in the maximization routine (if needed by method). } 
      \item{\code{maxit}}{ Maximum number of iterations to be used in the maximization routine (if needed by method). } 
      \item{\code{abstol}}{ Optional convergence tolerance for the maximization routine.  } 
      \item{\code{iter.V0}}{ The value of V0 to be used in place of 0 when x0 is treated as fixed and V0=0. See manual for discussion of initial state variance. }
      \item{\code{trace}}{ A positive integer.  If not 0, a record will be created during maximization iterations (what's recorded depends on method of maximization). }
      \item{\code{MCInit}}{ If TRUE, do a Monte Carlo search of the initial condition space.} 
      \item{\code{numInits}}{ Number of random initial value draws if MCInit=TRUE (ignored otherwise). } 
      \item{\code{numInitSteps}}{ Number of EM iterations for each random initial value draw if MCInit=TRUE (ignored otherwise). }
      \item{\code{boundsInits}}{ Bounds on the uniform distributions from which initial values will be drawn if MCInit=TRUE (ignored otherwise). }
      \item{\code{silent}}{ Suppresses printing of progress bar, error messages and convergence information. } 
      \item{\code{silent}}{ Suppresses printing of progress bar, error messages and convergence information. }       
    }
  }
  \item{\code{method}}{ A string specifying the estimation method. MARSS 1.0 allows "kem" for Kalman-EM and "BFGS" for quasi-Newton. }
  \item{\code{par}}{ A list with 8 matrices of estimated parameter values Z, A, R, B, U, Q, x0, V0. }
  \item{\code{kf}}{ A list containing Kalman filter/smoother output.}
  \item{\code{numIter}}{ Number of iterations which were required for convergence. }
  \item{\code{convergence}}{ Convergence status and errors. 0 means converged successfully.  Anything else means an error or warning.}
  \item{\code{logLik}}{ Log-likelihood. }
}
}

\value{
  TRUE if no problems; otherwise a message describing the problems.
}

\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{marssm}} \code{\link{MARSSkem}} 
}


