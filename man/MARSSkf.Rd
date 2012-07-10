\name{MARSSkf}
\alias{MARSSkf}
\alias{MARSSkfas}
\title{ Kalman Filtering and Smoothing }
\description{
  Implements the Kalman filter/smoother for MARSS models. This is a base function in the \code{\link{MARSS-package}}.  \code{MARSSkf()} is a native R implementation that allows x0 to be at t=0 or t=1 (data starts at t=1).  \code{MARSSkfas()} uses the Kalman filter-smoother functions in the KFAS package.  This requires that x0 be at t=1 (kf.x0="x10"), but the functions use more efficient and more stable algorithms.  Exact diffuse priors are also allowed in the KFAS Kalman filter function.
}
\usage{
MARSSkf( MLEobj )
}
\arguments{
  \item{ MLEobj }{ A \code{\link{marssMLE}} object with the \code{par} element of estimated parameters, \code{model} element with the model description and data, and \code{control} element for the fitting algorithm specificitions.  \code{control$debugkf} specifies that detailed error reporting will be returned.  \code{model$diffuse=TRUE} specifies that a diffuse prior be used (for \code{\link{MARSSkfas}} only. See KFAS documentation. V0 should be non-zero since the diffuse prior variance is V0*kappa, where kappa goes to infinity.}
  }
\details{
For state space models, the Kalman filter and smoother provide optimal (minimum mean square error) estimates of the hidden states. The Kalman filter is a forward recursive algorithm which computes estimates of the states x(t) conditioned on the data up to time t. The Kalman smoother is a backward recursive algorithm which starts at time T and works backwards to t = 1 to provide estimates of the states conditioned on all data.
 
The data may contain missing values.  In this case, the missing values are specified by \code{MLEobj$model$miss.value}. 

The expected value of the initial state, x0, is an estimated parameter (or treated as a prior).  This E(initial state) can be treated in two different ways.  One can treat it as x00, meaning E(x at t=0 | y at t=0), and then compute x10, meaning E(x at t=1 | y at t=0), from x00.    Or one can simply treat the initial state as x10, meaning E(x at t=1 | y at t=0).  The approaches lead to the same parameter estimates, but the likelihood is written slightly differently in each case and you need your likelihood calculation to correspond to how the initial state is treated in your model (either x00 or x10).  The EM algorithm in the MARSS package (\code{\link{MARSSkem}}) follows Shumway and Stoffer's derivation and uses x00, while Ghahramani et al uses x10.  The \code{MLEobj$model$tinitx} argument specifies whether the initial states (specified with x0 and V0) is at t=0 (\code{tinitx=0}) or t=1 (\code{tinitx=1}).

The Kalman filter/smoother code has been altered to allow the user to specify partially deterministic models (some or all elements of the Q diagonal equal to 0) and partially perfect observation models (some or all elements of the R diagonal equal to 0).  In addition, the code includes numerous checks to alert the user if matrices are becoming ill-conditioned and the algorithm unstable.  In addition, the likelihood computation returns the exact likelihood when there are missing values (rather than the approximate likelihood typically used).

All parameters may be time-varying.  See \code{\link{marssm}} for a description of marssm model objects and of how time-varying parameters are specified.
}
\value{
A list with the following components (m is the number of state processes). "V" elements are called "P" in Shumway and Stoffer.  
  \item{xtT}{ State first moment E[x(t) | y(1:T)] (m x T matrix). Kalman smoother output.}
  \item{VtT}{ State second moments E[x(t)x(t)'| y(1:T)] (m x m x T array). Kalman smoother output. P_t^T in S&S.}
  \item{Vtt1T}{ State lag-one second moments E[x(t)x(t-1)' | y(1:T)] (m x m x T). Kalman smoother output. P_{t,t-1}^T in S&S.}
  \item{x0T}{ Initial state estimate E[x(i) | y(1:T)] (m x 1). If control$kf.x0="x00", i=0; if ="x10", i=1. Kalman smoother output.}
  \item{V0T}{ Estimate of initial state covariance matrix E[x(i)x(i)' | y(1:T)] (m x m). If control$kf.x0="x00", i=0; if ="x10", i=1.  Kalman smoother output. P_0^T in S&S.}
  \item{J}{ (m x m x T) Kalman smoother output.  }
  \item{J0}{ J at init time (t=0 or t=1) (m x m x T). Kalman smoother output.  }
  \item{xtt}{ E[x(t) | y(1:t)] (m x T).  Kalman filter output.}
  \item{xtt1}{ E[x(t) | y(1:t-1)] (m x T).  Kalman filter output.}
  \item{Vtt}{ State second moment estimates, E[x(t)x(t)'| y(1:t)] (n x n x T). Kalman filter output. P_t^t in S&S.}
  \item{Vtt1}{ State second moment estimates E[x(t)x(t)' | y(1:t-1)] (m x m x T). Kalman filter output. P_t^{t-1} in S&S.}
  \item{Kt}{ Kalman gain (m x m x T). Kalman filter output. }
  \item{Innov}{ Innovations y(t) - E[y(t) | Y(t-1)] (n x T). Kalman filter output. }
  \item{Sigma}{ Innovations variances. Kalman filter output. }
  \item{logLik}{ Log-likelihood computed from \code{mssm.params} and innovations. }
  \item{errors}{ Any error messages. }
}
\references{
A. C. Harvey (1989).  Chapter 5, Forecasting, Structural Time Series Models and the Kalman Filter.  Cambridge University Press.

R. H. Shumway and D. S. Stoffer (2006).  Chapter 6, Time Series Analysis and its Applications.  Springer-Verlag, New York.

Ghahramani, Z. and Hinton, G.E. (1996) Parameter estimation for linear dynamical systems.  University of Toronto Technical Report CRG-TR-96-2.

The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.     

}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}} \code{\link{marssm}} \code{\link{MARSSkem}}
}

\examples{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  #you can use MARSS to construct a MLEobj
  #MARSS calls MARSSinits to construct default initial values
  MLEobj = MARSS(dat, fit=FALSE)
  #MARSSkf needs a marss MLE object with the par element set
  MLEobj$par=MLEobj$start
  #Compute the kf output at the params used for the inits 
  kfList = MARSSkf( MLEobj )
}

