\name{print.marssMLE}
\alias{print}
\alias{print.MARSS}
\alias{print.marssMLE}
\title{ Printing function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{print(marssMLE)}., where marssMLE is one's output from a \code{\link{MARSS}} call, will print out information on the fit.  However, \code{print} can be used to print a variety of information from a marssMLE object using the \code{what} argument in the print call.
}
\usage{
\method{print}{marssMLE}(x, digits = max(3, getOption("digits")-4), ..., what="fit", silent=FALSE)
}
\arguments{
  \item{x}{ A marssMLE object.  }
  \item{digits}{ Number of digits for printing.  }
  \item{...}{ Other arguments for print. }
  \item{what}{ What to print.  Default is "fit". 
  \itemize{
    \item{ "model"}{ The model parameters with names for the estimted parameters.  The output is customized by the form of the model that was fit.  This info is in \code{marssMLE$call}. }
    \item{ "par" }{ A vector of all the estimated values. }
    \item{ "xtT" or "states" }{ The estimated states. }
    \item{ "data" }{ The data. }
    \item{ "ytT" }{ The expected value of the data.  Returns the data if present and the expected value if missing. }
    \item{ "states.se" }{ The states standard errors. }
    \item{ "states.cis" }{ Approximate confidence intervals for the states. }
    \item{ parameter name }{ Returns the parameter matrix for that parameter with fixed values at their fixed values and the estimated values at their estimated values. }
  } }
  \item{silent}{ If TRUE, do not print just return the invisible object if print call is assigned.  See example.}
}

\value{
  A print out of information.  If you assign the print call to a value, then you can reference the output.  See the examples.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  
  print(MLEobj)
  
  print(MLEobj,what="model")
  
  print(MLEobj,what="par")
  
  cis=print(MLEobj,what="states.cis",silent=TRUE)
  cis$up95CI
  
  R=print(MLEobj, what="R")
  R
}