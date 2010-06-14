\name{MARSSapplynames}
\alias{MARSSapplynames}
\title{ Names for marssMLE Object Components }
\description{
  Puts names on the matrix components of \code{\link{marssMLE}} and \code{\link{marssm}} objects.     This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSapplynames(obj, Y.names = NA, X.names = NA, x0.names = NA, 
  V0.names = NA, U.names = NA, A.names = NA, R.names = NA, 
  Q.names = NA, B.names = NA, Z.names = NA,
  rows = TRUE, cols = TRUE)
}
\arguments{
  \item{obj}{ An object of class \code{\link{marssMLE}} or \code{\link{marssm}}. }
  \item{Y.names}{ Vector of names for observed time series. }
  \item{X.names}{ Vector of names for the hidden state trajectories. }
  \item{x0.names}{ Vector of names for \code{MLEobj$par$x0}, \code{MLEobj$start$x0}, 
\code{MLEobj$model$fixed$x0} and \code{MLEobj$model$free$x0}. }
  \item{V0.names}{ Names for \code{MLEobj$par$V0} }
  \item{U.names}{ Names for \code{MLEobj$par$U}  }
  \item{A.names}{ Names for \code{MLEobj$par$A}  }
  \item{R.names}{ Names for \code{MLEobj$par$R} }
  \item{Q.names}{ Names for \code{MLEobj$par$Q}  }
  \item{B.names}{ Names for \code{MLEobj$par$B}  }
  \item{Z.names}{ Row and col names for \code{MLEobj$par$Z}}
  \item{rows}{ Add row names? }
  \item{cols}{ Add column names to B, Q, R and V0 matrices? }
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
\examples{
  dat = t(harborSealWA)
  dat = dat[2:nrow(dat),]
  kem = MARSS(dat)
  kemNamed = MARSSapplynames(kem)
}

