\name{popWrap}
\alias{popWrap}
\title{ Wrapper Objects }
\description{
  Wrapper objects are used by the function \code{\link{MARSS}} in the package \code{\link{MARSS-package}}.  \code{popWrap} creates a wrapper object containing specifications and options for estimation 
  of a multivariate autoregressive state-space model.
}
\usage{
popWrap(y, allowed,
    inits=NULL,
    model=NULL,
    fixed=NULL, free=NULL, 
    miss.value=NULL,
    control=NULL,
    method=NULL,   
    silent=FALSE)
}
\arguments{
  \item{y}{ A matrix (not dataframe), observations (rows) x time steps (columns), of data. }
  \item{allowed}{ Allowed model structures. This is determined by the method used for parameter estimation.  \code{allowed$kem} (specified in \code{MARSSsettings} is the set of allowed model structures for the Kalman-EM algorithm under MARSS 1.0. \code{allowed$BFGS} are the allowable model structures for BFGS.}
  \item{inits}{ List with up to 8 matrices Z, A, R, B, U, Q, x0, V0, specifying initial values for parameters to use in iterative ML estimation algorithms, such as Kalman-EM and quasi-Newton methods. These are ignored if the specified parameter is not being estimated.  If not passed in by the user, \code{\link{MARSS}} creates generic inits using \code{\link{MARSSinits}}.
    \describe{
    \item{\code{B}}{ Initial value(s) for B parameter (length 1 or m x m). If length 1, \code{inits$B} constructed as \code{diag(value, m)}.}
    \item{\code{U}}{ Initial value(s) for U parameter (length 1 or m x 1). If length 1, \code{inits$U} constructed as \code{matrix(value, nrow=m, ncol=1)}.}
    \item{\code{Q}}{ Initial value(s) for Q parameter (length 1 or m x m). If length 1, \code{inits$Q} constructed as \code{diag(value, m)}. }
    \item{\code{Z}}{ Initial value(s) for Z parameter (n x m). Ignored in MARSS 1.0; included for MARSS 2.0.}
    \item{\code{A}}{ Initial value(s) for A parameter (length 1 or n x 1). If length 1, \code{inits$U} constructed as \code{matrix(value, nrow=n, ncol=1)}.}
    \item{\code{R}}{ Initial value(s) for R parameter (length 1 or n x n). If length 1, \code{inits$R} constructed as \code{diag(value, n)}.}
    \item{\code{x0}}{ Initial value(s) for x0 parameter (length 1 or m x 1). If length 1, \code{inits$x0} constructed as \code{matrix(value, nrow=m, ncol=1)}.}
    \item{\code{V0}}{ Initial variance(s) for hidden states (length 1 or m x m). Ignored in MARSS 1.0; included for forward compatibility.}
    }
    }
  \item{model}{ Model specification as a list. See \code{\link{MARSS}} for details. }
  \item{fixed}{ Optional model specification using matrices of fixed and free parameters. See user guide for details. }
  \item{free}{ Optional model specification using matrices of fixed and free parameters. See user guide for details. }
  \item{miss.value}{ How are missing values represented in the data? }
  \item{method}{ The method used for estimation.  This is needed for setting default values for \code{control}. }
  \item{control}{ Control options for maximization algorithms. 
    \describe{
      \item{\code{minit}}{ Minimum number of EM iterations. } 
      \item{\code{maxit}}{ Maximum number of EM iterations. } 
      \item{\code{abstol}}{ Optional tolerance for log-likelihood change.  If log-likelihood changes less than this amount relative to the previous iteration, the algorithm exits. } 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{ Positive integer. If not 0, a record will be created of each variable over the maximization iterations. The information recorded depends on the maximization method.}
      \item{\code{safe}}{ If TRUE, \code{\link{MARSSkem}} will rerun \code{\link{MARSSkf}} after each individual parameter update rather than only after all parameters are updated.  }
      \item{\code{MCInit}}{ Use Monte Carlo initialization? See discussion in \code{\link{MARSSkem}} and \code{\link{MARSSmcinit}}. } 
      \item{\code{numInits}}{ Number of random initial value draws. } 
      \item{\code{numInitSteps}}{ Number of EM iterations for each random initial value draw. }
      \item{\code{boundsInits}}{ Bounds on the uniform distributions from which initial values will be drawn. Note that bounds for the covariance matrices Q and R, which require positive values, are specified in logs. }    
    }
  }
  \item{silent}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }
}
\details{
Wrapper functions e.g. \code{\link{MARSS}} call \code{popWrap()} to create a 'popWrap' object, then \code{\link{is.marssm}} to coerce this object to class 'marssm' for the estimation function.  The popWrap() function calls \code{\link{checkPopWrap}} to check user inputs.

If arguments \code{inits}, \code{model}, or \code{control} are not provided by the user, they will be set by the \code{alldefaults[[method]]} object specified in \code{MARSSsettings}.  Argument \code{model} specifies the model structure using a list of matrices; see \code{\link{MARSS}} or the user guide for instructions on how to specify model structure.
}
\value{
  An object of class 'popWrap'.
  \item{data}{ Data supplied by user. }
  \item{m}{ Number of hidden state trajectories. }
  \item{model}{ A list with up to 8 elements Z, A, R, B, U, Q, x0, V0 (unless some of these are specified in "fixed"). See \code{\link{MARSS}} for details on what values are allowed. }
  \item{fixed}{ A list with (up to) 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{free}{ A list with (up to) 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{inits}{ A list specifying initial values for parameters to be used at iteration 1 in iterative maximum-likelihood algorithms. }
  \item{miss.value}{ Specifies missing value representation. }
  \item{method}{ The method used for estimation. }
  \item{control}{ See Arguments. }
}

\author{ 
  Kellie Wills, NOAA, Seattle, USA.

  kellie(dot)wills(at)noaa(dot)gov   
}
\seealso{ 
  \code{\link{MARSS}} \code{\link{marssm}} \code{\link{checkPopWrap}} \code{\link{as.marssm}} 
}
\examples{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  wrapperObj = popWrap(dat, allowed=MARSS:::allowed$kem, method="kem")
}

