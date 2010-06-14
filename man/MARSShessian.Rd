\name{MARSShessian}
\alias{MARSShessian}
\title{ MARSS Parameter Variance-Covariance Matrix from the Hessian Matrix}
\description{
  Calculates approximate parameter variance-covariance matrix and appends it to a \code{\link{marssMLE}} object.  This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSShessian(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
}
\details{
  Uses \code{\link{fdHess}} from package \code{\link{nlme}} to numerically estimate the Hessian matrix
   (the matrix of partial 2nd derivatives of the parameter estimates).  Hessian CIs are based on 
   the asymptotic normality of ML estimates under a large-sample approximation. 
}
\value{
  \code{MARSShessian()} returns the \code{\link{marssMLE}} object passed in along with additional components
   \code{Hessian, gradient, parMean} and \code{parSigma} computed by the \code{MARSShessian} function.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSparamCIs}} \code{\link{marssMLE}} }
\examples{
  \dontrun{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  MLEobj = MARSS(dat)
  MLEobj.w.hessian = MARSShessian(MLEobj)
  }
}

