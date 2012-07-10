\name{MARSSvectorizeparam}
\alias{MARSSvectorizeparam}
\title{ Vector to Parameter Matrix Conversion }
\description{
  Converts \code{MLEobj$par} to a vector of the estimated parameter elements and vice versa.  This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSvectorizeparam(MLEobj, parvec = NA)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
  \item{parvec}{ NA or a vector. See Value. }
}
\details{
  Utility function to generate parameter vectors for optimization functions, and to set \code{MLEobj$par} using a vector of parameter values (only the estimated values). 
}
\value{
  If parvec=NA, a vector of estimated parameters. Otherwise, a \code{\link{marssMLE}} object with \code{$par} set by parvec.
}
\author{ 
  Eli Holmes and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ \code{\link{marssMLE}} }
\examples{
  dat = t(harborSealWA)
  dat = dat[2:4,]
  kem = MARSS(dat)
  paramvec = MARSSvectorizeparam(kem)
  paramvec
}

