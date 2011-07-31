\name{MARSShatyt}
\alias{MARSShatyt}
\title{ Compute Expected Value of Y,YY, and YX}
\description{
  Computes the expected value of random variables involving Y for the EM algorithm. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSShatyt(y, parList, kfList, missing.matrix = NULL, miss.value= NULL)
}
\arguments{
  \item{y}{ A matrix (not dataframe), sites (rows) x years (columns). See Details regarding handling of missing values. }
  \item{parList}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0, specifying parameter values. An example is the \code{par} element in a \code{\link{marssMLE}} object.}
  \item{kfList}{ A list with the output from \code{\link{MARSSkf}}.}
  \item{missing.matrix}{ Optional matrix specifying which observations are missing. See Details. }
  \item{miss.value}{ How are missing values represented in the data? Either miss.value or missing.matrix must be supplied. If both are supplied, then miss.value will be ignored with no warning(!). }
  }
\details{
For state space models, \code{MARSShatyt()} computes the expectations involving Y.  If Y is completely observed, this entails simply replacing Y with the observed y.  When Y is only partially observed, the expectation involves the conditional expectation of a multivariate normal. 
 
Missing values in the data may be handled in one of two ways:
1. Missing values may be replaced with zeroes prior to passing to \code{MARSSkf()}. Argument \code{missing.matrix} must then be a matrix of the same dimensions as the data, with 0 in the positions of observed values and 1 in the positions of missing values. 
2. Data containing missing values may be passed in. Argument \code{miss.value} must then be the code used to represent missing values. 
The function requires that you specify either a missing matrix or a miss.value.  If there are no missing values, just set \code{miss.value} to a value that is not in your data (like NA or -99).

}
\value{ 
A list with the following components (n is the number of state processes). Names ending in "T" are estimates from the Kalman smoother; J is also smoother output. Other components are output from the Kalman filter.   
  \item{ytT}{ Estimates E[Y(t) | Y(1)=y(1)] (n x T matrix). }
  \item{OtT}{ Estimates E[Y(t)t(Y(t) | Y(1)=y(1)] (n x n x T array). }
  \item{yxtT}{ Estimates E[Y(t)t(X(t) | Y(1)=y(1)] (n x m x T array). }
  \item{errors}{ Any error messages due to ill-conditioned matrices. }
  \item{ok}{ (T/F) Whether errors were generated. }
}
\references{
Holmes, E. E. (2010) Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models. Technical report. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}} \code{\link{marssm}} \code{\link{MARSSkem}}
}


