\name{residuals.marssMLE}
\alias{residuals}
\alias{residuals.marssMLE}
\alias{residuals.MARSS}
\alias{residuals.marss}
\title{ MARSS standardized residuals }
\description{
  Calculates the standardized (or auxilliary) residuals sensu Harvey, Koopman and Penzer (1998).  This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
\method{residuals}{marssMLE}(object,...)
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}. 
  This object must have \code{$par}, \code{$model} and \code{$kf} elements containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. This may require setting \code{control$trace=1} in the \code{MARSS} call so that the \code{kf} list is attached to the MLE object.}
    \item{...}{ Not used.}
}
\details{
  Uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998) to compute the standardized model residuals.  
}
\value{
A list with the following components  
  \item{vthat}{ The smoothed model residuals y(t)-E(y(t)|y(1:T),Theta), where Theta is the set of model parameters (Z, A, R, etc).  Sometimes called the smoothations. This is different than the Kalman filter innovations which are y(t)-E(y(t)|y(1:t-1),Theta).}
  \item{wthat}{ The smoothed stated residuals E(x(t)|y(1:T))-E(x(t)|E(x(t-1)|y(1:T))).}
  \item{et}{ The model residuals as a (n+m) x TT matrix with vthat on top and wthat below.    et is hat(eta_t) on page 112 of Harvey, Koopman and Penzer (1998).  }
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
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  
  #not standardized model residuals
  residuals(MLEobj)$model.residuals

  #standardized (by variance) model & state residuals
  residuals(MLEobj)$std.residuals
}
\references{
Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified approach. Advances in Econometrics 13: 103-144  (see page 112).

Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166. (see pages 147-148).
}

