\name{find.degenerate}
\alias{find.degenerate}
\title{
Find degenerate variance parameters
}
\description{
A helper function to find degenerate variance parameters in Kalman-EM estimates in the package \code{\link{MARSS-package}}.
}
\usage{
find.degenerate(MLEobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MLEobj}{
  An object of class marssMLE as output by \code{\link{MARSSkem}}. Typically after a call to \code{\link{MARSS}}: \code{MARSS(data,method="kem")}.
}
}
\details{
This function plots the log of the absolute value of the variance element (on the diagonal) against the log iteration number.  Such log-log plots are
commonly used to assess convergence in iterative routines.   In state-space models with both process and non-process (observation) variance, it is
entirely possible that the highest likelihood occurs when one of the variance element is zero.  In this case, that element is degenerate (=0).  Since
the likelihood computations will generate an error, this means you will not be able to compute the true likelihood.  Technically, it is certaintly
possible to compute the likelihood, but the code in \code{\link{MARSSkf}} will throw an error if you try to set one of the variances to zero. 
}
\value{
  A plot where each diagonal variance element appears in a panel.  Converged elements with have a flat log-log plot.  Degenerate variance elements will
  have a declining log-log plot.  Typically, degenerate variance parameters will show a slope of ca -1.
}
\author{
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{
\code{ \link{MARSSkem}} }
\examples{
  dat = t(harborSealnomiss)
  #shorten the data set to 10 yrs in order to produce a degenerate solution
  #in this case there is not enough data to estimate all the Q variances
  dat = dat[2:nrow(dat),1:10]
  #force 200 interations to get a nice log(param):log(iteration) plot
  MLEobj = MARSS(dat, control=list(minit=200))
  find.degenerate(MLEobj)
}
