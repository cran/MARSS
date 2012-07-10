\name{MARSSapplynames}
\alias{MARSSapplynames}
\title{ Names for marssMLE Object Components }
\description{
  Puts names on the matrix components of \code{\link{marssMLE}} and \code{\link{marssm}} objects.  This is a utility function in the \code{\link{MARSS-package}} and is not directly accessible.  Use \code{MARSS:::MARSSapplynames()} to access.
}
\usage{
MARSSapplynames(obj, Y.names = NA, X.names = NA)
}
\arguments{
  \item{obj}{ An object of class \code{\link{marssMLE}} or \code{\link{marssm}}. }
  \item{Y.names}{ Vector of names for observed time series. }
  \item{X.names}{ Vector of names for the hidden state trajectories. }
}
\details{
Default behavior will use names rownames of data (if available, and if not Y1, Y2, ...) for the Ys and all matrices that are n x ... and use X1, X2, ... for the Xs and all matrices that are m x ... .
}
\value{
  The object passed in, with row and column names on matrices as specified.
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, kellie(dot)wills(at)noaa(dot)gov
}
\seealso{ \code{\link{marssMLE}} \code{\link{marssm}} }


