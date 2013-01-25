\name{coef.marssMLE}
\alias{coef}
\alias{coef.MARSS}
\alias{coef.marssMLE}
\alias{coef.marss}
\title{ Coefficient function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{coef(marssMLE)}., where marssMLE is one's output from a \code{\link{MARSS}} call, will print out the estimated parameters.  The default output is a list of the estimated parameters, however \code{coef} can be altered using the \code{what} argument in the function call.
}
\usage{
\method{coef}{marssMLE}(object, ..., type="list", form=NULL)
}
\arguments{
  \item{object}{ A marssMLE object.  }
  \item{...}{ Other arguments for coef. }
  \item{type}{ What to print.  Default is "list".  If you input type as a vector, coef returns a list of output. See examples.
  \itemize{
    \item{ "list" }{ A list of only the estimated values in each matrix. Each model matrix has it's own list element.}
    \item{ "vector" }{ A vector of all the estimated values in each matrix. }
    \item{ "matrix" }{ Returns the parameter matrix for that parameter with fixed values at their fixed values and the estimated values at their estimated values. }
  } }
  \item{form}{ By default, coef uses the model form specified in the call to \code{\link{MARSS}} to determine how to display the coefficients.  This information is in \code{ marssMLE$call$form }, however you can specify a different form.  \code{form="marss"} should always work since this is the model form in which the model objects are stored (in \code{marssMLE$model}).}
}

\value{
  A list of the estimated parameters for each model matrix.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{ 
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  
  coef(MLEobj)
  coef(MLEobj,type="vector")
  coef(MLEobj,type="matrix")
}