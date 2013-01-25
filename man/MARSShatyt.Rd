\name{MARSShatyt}
\alias{MARSShatyt}
\title{ Compute Expected Value of Y,YY, and YX}
\description{
  Computes the expected value of random variables involving Y for the EM algorithm. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSShatyt( MLEobj )
}
\arguments{
  \item{ MLEobj }{ A \code{\link{marssMLE}} object with the \code{par} element of estimated parameters, \code{model} element with the model description and data. }
  }
\details{
For state space models, \code{MARSShatyt()} computes the expectations involving Y.  If Y is completely observed, this entails simply replacing Y with the observed y.  When Y is only partially observed, the expectation involves the conditional expectation of a multivariate normal. 
 
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
Holmes, E. E. (2012) Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models. Technical report. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112. Type \code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}} \code{\link{marssm}} \code{\link{MARSSkem}}
}


