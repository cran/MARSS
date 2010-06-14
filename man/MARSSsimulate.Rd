\name{MARSSsimulate}
\alias{MARSSsimulate}
\title{ Simulate Data from a MARSS Model and Parameter Estimates }
\description{
  Generates simulated data from a MARSS model with specified parameter estimates. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSsimulate(parList, tSteps = 100, nsim = 1, silent = TRUE, 
  miss.loc = NULL, miss.value = NULL)
}
\arguments{
  \item{parList}{ A list of parameter matrices structured like the \code{$par} element of an object of class \code{\link{marssMLE}}. }
  \item{tSteps}{ Number of time steps in each simulation. }
  \item{nsim}{ Number of simulated data sets to generate. }
  \item{silent}{ Suppresses progress bar. }
  \item{miss.loc}{ Optional matrix specifying where to put missing values. See Details. }
  \item{miss.value}{ Code representing missing values in \code{miss.matrix}. See Details. }
}
\details{
  Argument \code{miss.loc} is an array of dimensions n x tSteps x nsim, specifying where to put missing values
  in the simulated data. Locations where missing data appear are specified using the \code{miss.value}; non-missing values can be specified by any other numeric value.  If the locations of the missing values are the same for all simulations, \code{miss.loc} can be a matrix of dim=c(n, tSteps) (the original data for example). If \code{miss.loc} is passed in, \code{miss.value} must be specified.  The default is that there are no missing values.  If one's data has missing values in it and one want to replicate those locations in the simulated data, \code{miss.loc} can simply be set to the original data (see examples).
}
\value{
  \item{sim.states}{ Array (dim m x tSteps x nsim) of state processes simulated from parameter estimates. }
  \item{sim.data}{ Array (dim n x tSteps x nsim) of data simulated from parameter estimates.}
  \item{par}{ The list of parameter matrices from which the data were simulated. }
  \item{miss.loc}{ Matrix identifying where missing values are located. }
  \item{tSteps}{ Number of time steps in each simulation. }
  \item{nsim}{ Number of simulated data sets generated. }
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{marssm}} \code{\link{marssMLE}} \code{\link{MARSSboot}} 
}
\examples{
  #do a parametric bootstrap.  
  #Same length as original data and same location of missing data
  d = harborSeal[,2:ncol(harborSeal)]
  dat = t(d)
  MLEobj = MARSS(dat)
  sim.obj = MARSSsimulate(parList=MLEobj$par, tSteps=dim(d)[1], nsim=10)
}

