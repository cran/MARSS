\name{MARSSkemcheck}
\alias{MARSSkemcheck}
\title{ Model Checking for MLE objects passed to MARSSkem }
\description{
  This is a helper function in the \code{\link{MARSS-package}} that checks that the model can be handled by the \code{\link{MARSSkem}} algorithm.
}
\usage{
MARSSkemcheck(modelObj)
}
\arguments{
  \item{modelObj}{ An object of class \code{\link{marssm}}. }
}
\value{
  A list with of the model elements A, B, Q, R, U, x0, Z, V0 specifying the structure of the model using text strings).
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{marssm}} \code{\link{MARSSkem}} }

