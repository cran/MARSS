\name{MARSSsimulate}
\alias{MARSSsimulate}
\title{ Simulate Data from a MARSS Model and Parameter Estimates }
\description{
  Generates simulated data from a MARSS model with specified parameter estimates. This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSsimulate(MLEobj, tSteps = 100, nsim = 1, silent = TRUE, 
  miss.loc = NULL)
}
\arguments{
  \item{MLEobj}{ A marss MLE model object \code{\link{marssMLE}}. }
  \item{tSteps}{ Number of time steps in each simulation. }
  \item{nsim}{ Number of simulated data sets to generate. }
  \item{silent}{ Suppresses progress bar. }
  \item{miss.loc}{ Optional matrix specifying where to put missing values. See Details. }
}
\details{
  Argument \code{miss.loc} is an array of dimensions n x tSteps x nsim, specifying where to put missing values
  in the simulated data. If missing, this would be constructed using \code{modelObj$data} and \code{modelObj$miss.value}.  If the locations of the missing values are the same for all simulations, \code{miss.loc} can be a matrix of dim=c(n, tSteps) (the original data for example). If \code{miss.loc} is passed in, missing values must be specified as \code{modelObj$miss.value}.  The default, if \code{miss.loc} is left off, is that there are no missing values even if \code{modelObj$data} has missing values. 
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
  d = harborSeal[,c(2,11)]
  dat = t(d)
  MLEobj = MARSS(dat)
  sim.obj = MARSSsimulate(MLEobj, tSteps=dim(d)[1], nsim=10)
}

