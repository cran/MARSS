\name{MARSSparamCIs}
\alias{MARSSparamCIs}
\title{ Confidence Intervals for MARSS Parameters }
\description{
  Computes confidence intervals for the maximum-likelihood estimates of MARSS model parameters.  This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSparamCIs(MLEobj, method = "hessian", alpha = 0.05, nboot=1000)
}

\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. Must have a \code{$par} element containing the MLE parameter estimates. }
  \item{method}{ Method for calculating the standard errors: "hessian", "parametric", and "innovations" implemented currently. }
  \item{alpha}{ alpha level for the 1-alpha confidence intervals. }
  \item{nboot}{ Number of bootstraps to use for "parametric" and "innovations" methods. }
}
\details{
  Approximate confidence intervals (CIs) on the model parameters may be calculated from the Hessian matrix (the matrix of partial 2nd derivatives of the parameter estimates) or parametric or non-parametric (innovations) bootstrapping using \code{nboot} bootstraps.  The Hessian CIs are based on the asymptotic normality of MLE parameters under a large-sample approximation.  Bootstrap estimates of parameter bias are reported if method "parametric" or "innovations" is specified.
}
\value{
  \code{MARSSparamCIs} returns the \code{\link{marssMLE}} object passed in, with additional components \code{par.se, par.upCI, par.lowCI, par.CI.alpha, par.CI.method, par.CI.nboot} and \code{par.bias} (if method is "parametric" or "innovations").
}
\references{ 
Holmes, E. E. and E. J. Ward. 2010. Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112; this is the user manual accesses via \code{show.doc(MARSS, manual)}
 }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{MARSSboot}} \code{\link{MARSSinnovationsboot}} \code{\link{MARSShessian}}
 }
\examples{
  dat = t(harborSealWA)
  dat = dat[2:4,]
  kem = MARSS(dat, constraint=list(Z=factor(c(1,1,1)), 
     R="diagonal and unequal"))
  kem.with.CIs.from.hessian = MARSSparamCIs(kem)
  kem.with.CIs.from.hessian
}

