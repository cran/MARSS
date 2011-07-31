\name{MARSSresids}
\alias{MARSSresids}
\title{ MARSS standardized residuals }
\description{
  Calculates the standardized residuals sensu Harvey, Koopman and Penzer (1998).  This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSresids(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have \code{$par}, \code{$model} and \code{$kf} elements containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
}
\details{
  Uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998) to compute the standardized model residuals.  
}
\value{
A list with the following components  
  \item{et}{ The model residuals as a (n+m) x TT matrix with v_t on top and w_t below.  This is hat(eta_t) on page 112 of Harvey, Koopman and Penzer (1998).  }
  \item{var.et}{ The variance of the model residuals as a (n+m) x (n+m) x TT matrix. This is var(hat(eta_t)).}
  \item{std.et}{ The standardized model residuals as a (n+m) x TT matrix. This is \code{et} divided by the square root of \code{var.et} --- although the code is using the matrix equivalent of that equation.}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSkem}} \code{\link{marssMLE}} }
\examples{
  dat = t(harborSeal)
  dat = dat[2:4,]
  MLEobj = MARSS(dat)
  std.resids = MARSSresids(MLEobj)$std.et
}
\references{
Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified approach. Advances in Econometrics 13: 103-144  (see page 112).

Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166. (see pages 147-148).
}

