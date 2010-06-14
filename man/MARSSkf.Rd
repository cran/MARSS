\name{MARSSkf}
\alias{MARSSkf}
\title{ Kalman Filtering and Smoothing }
\description{
  Implements the Kalman filter/smoother for MARSS models. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSkf(y, parList, missing.matrix = "not given", 
  miss.value = "not given", init.state="x10", debugkf=FALSE)
}
\arguments{
  \item{y}{ A matrix (not dataframe), sites (rows) x years (columns). See Details regarding handling of missing values. }
  \item{parList}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0, specifying parameter values. An example is the \code{par} element in a \code{\link{marssMLE}} object.}
  \item{missing.matrix}{ Optional matrix specifying which observations are missing. See Details. }
  \item{miss.value}{ How are missing values represented in the data? Either miss.value or missing.matrix must be supplied. If both are supplied, then miss.value will be ignored with no warning(!). }
  \item{init.state}{ Is the initial state, x0, treated as E(x) at time t=0 (\code{init.state="x00"}) or E(x) at t=1 (\code{init.state="x10"})? Default is \code{init.state="x10"} (note! for Shumway and Stoffer's Kalman filter use \code{init.state="x00"}).  See Details. }
  \item{debugkf}{ Return detailed error messages?}
  }
\details{
For state space models, the Kalman filter and smoother provide optimal (minimum mean square error) estimates of the hidden states. The Kalman filter is a forward recursive algorithm which computes estimates of the states x(t) conditioned on the data up to time t. The Kalman smoother is a backward recursive algorithm which starts at time T and works backwards to t = 1 to provide estimates of the states conditioned on all data.
 
Missing values in the data may be handled in one of two ways:
1. Missing values may be replaced with zeroes prior to passing to \code{MARSSkf()}. Argument \code{missing.matrix} must then be a matrix of the same dimensions as the data, with 0 in the positions of observed values and 1 in the positions of missing values. 
2. Data containing missing values may be passed in. Argument \code{miss.value} must then be the code used to represent missing values. 
The function requires that you specify either a missing matrix or a miss.value.  If there are no missing values, just set \code{miss.value} to a value that is not in your data (like -99).

The expected value of the initial state, x0, is an estimated parameter (or treated as a prior).  This E(initial state) can either be treated in two different ways.  One can treat it as x00, meaning E(x at t=0 | y at t=0), and then compute x10, meaning E(x at t=1 | y at t=0), from x00.    Or one can simply treat the initial state as x10, meaning E(x at t=1 | y at t=0).  The approaches are equivalent, but the likelihood is written slightly differently in each case (sum over 2 to T versus 1 to T) and you need your likelihood calculation to correspond to how the initial state is treated in your model (either x00 or x10).  The EM algorithm in the MARSS package (\code{\link{MARSSkem}}) follows Ghahramani's derivation and uses x10, while Shumway and Stoffer uses x00.  The \code{init.state} argument specifies whether x00 (\code{init.state="x00"}) or x10 (\code{init.state="x10"}) is used.  The default is \code{init.state="x10"}.
}
\value{
A list with the following components (n is the number of state processes). Names ending in "T" are estimates from the Kalman smoother; J is also smoother output. Other components are output from the Kalman filter.   
  \item{xtT}{ State estimates E[x(t) | y(1:T)] (n x T matrix). }
  \item{VtT}{ State covariances E[V(t) | y(1:T)] (n x n x T array). }
  \item{Vtt1T}{ Conditional error covariances E[V(t,t-1) | y(1:T)] (n x n x T). }
  \item{x0T}{ Initial state mean estimates (n x 1). }
  \item{V0T}{ Estimate of initial state covariance matrix (n x n). }
  \item{J}{ n x n x T  }
  \item{xtt1}{ Forecasts E[x(t) | Y(t-1)] (n x T).  }
  \item{Vtt}{ State covariance estimates E[V(t) | y(1:t)] (n x n x T). }
  \item{Vtt1}{ Conditional error covariances E[V(t,t-1) | y(1:t)] (n x n x T). }
  \item{Kt}{ Kalman gain (n x n x T). }
  \item{Innov}{ Innovations y(t) - E[y(t) | Y(t-1)] (n x T) (n x n x T). }
  \item{Sigma}{ Innovations variances. }
  \item{logLik}{ Log-likelihood computed from \code{mssm.params} and innovations. }
  \item{errors}{ Any error messages due to ill-conditioned matrices. }
}
\references{
P. J. Brockwell and R. A. Davis (1991).  Time Series: Theory and Methods.
  
A. C. Harvey (1989).  Chapter 5, Forecasting, Structural Time Series Models and the Kalman Filter.  Cambridge University Press.

R. H. Shumway and D. S. Stoffer (2006).  Chapter 6, Time Series Analysis and its Applications.  Springer-Verlag, New York.

The user manual:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   type \code{show.doc(MARSS, manual)} to see.          

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
  #Compute the kf output at the params used for the inits 
  kfList = MARSSkf(dat, MLEobj$start, miss.value=-99)
}

