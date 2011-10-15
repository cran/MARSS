\name{MARSSLLprofile}
\alias{MARSSLLprofile}
\title{ Log-likelihood profiles for MARSS Parameters }
\description{
  Computes log-likelihood profiles for the maximum-likelihood estimates of MARSS model parameters.  This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSLLprofile(MLEobj, param=NULL, x=NULL, LLlim=3, pstep=0.01, 
               max.steps=20, plot=TRUE)
}

\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. Must have a \code{$par} element containing the MLE parameter estimates. }
  \item{param}{ A vector of parameter names.  Must match those output from \code{\link{MARSSvectorizeparam}} and output when a \code{\link{marssMLE}} object is printed. If you leave this off (leave NULL), then profiles will be computed for all free variables.}
  \item{x}{ An optional vector of parameters values at which to comput the log-likelihood. If x is not specified then LLlim must be. If you pass in \code{x}, then \code{LLlim} and \code{pstep} will be ignored. If you want a different \code{x} for each parameter, then pass in a list: \code{x=list(U=c(...), Q=c(...))}.  Any left off parameters will use \code{LLlim} and \code{pstep}.}
  \item{LLlim}{ If passed in, the LL profile will go from max(LL)-LLlim to max(LL) to max(LL)-LLlim. The range of parameter values will keep being stepped upward or downward to reach these limits.}
  \item{pstep}{ The value of step=parameter(x(1))-parameter(x(2)) to use when steping towards LLlim.  If the parameters are from Q, R or V0, then step is on the log-scale: step=log(parameter(x(1)))-log(parameter(x(2))). If you want different step values for different parameters, pass this in as a list: \code{pstep=list(U=0.01, Q=.1, R=.2)}. Any parameter names left off will use the default step of 0.01.}
  \item{max.steps}{ Since it is not guaranteed that LLlim will be reached, max.steps is used to specify the maximum number of steps to do.}
  \item{plot}{ Whether to produce a plot.  If FALSE, then just the x values and LL values are returned. }
  }
\details{
 Computes log-likelihood profiles for the free parameters specified in \code{param} for a fitted \code{\link{marssMLE}} object using the method \code{MLEobj$method}. A red line is plotted at max.LL-1.92 (corresponding to the 1 degree of freedom chi-square value) If you want to change the \code{control} values used for computing the log-likelihoods (to say speed things up), then change the \code{$control} element of \code{MLEobj}.  If you have a lot of free parameters, then just call \code{MARSSLLprofile} with a few parameters at a time.  Otherwise the function will try to plot all of them in one plot and you'll get a plot margin error.
}
\value{
  \code{MARSSLLprofile} returns a list of x-values and LL profile values for each parameter in \code{param}.
}
\references{ 
Holmes, E. E. and E. J. Ward. 2010. Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112; this is the user guide accessed via \code{RShowDoc("UserGuide",package="MARSS")}
 }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{MARSSparamCIs}}
 }
\examples{
  dat = t(harborSealnomiss)
  dat = dat[2:3,]
  kem = MARSS(dat, model=list(Z=factor(c(1,1)), 
     R="diagonal and equal",Q=matrix(0.01)))
  profLL = MARSSLLprofile(kem,c("R.1","U.1"),pstep=list(U=0.01,R=0.1))
}

