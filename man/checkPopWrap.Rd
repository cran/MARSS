\name{checkPopWrap}
\alias{checkPopWrap}
\title{ Check Arguments to popWrap() }
\description{
  Checks inputs for wrapper object (\code{\link{popWrap}}) creation to ensure that the wrapper object can be handled by \code{\link{as.marssm}}.  This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
checkPopWrap(wrapperObj, wrapper.el, allowed, silent=FALSE)
}
\arguments{
  \item{wrapperObj}{ An object of class \code{\link{popWrap}}. }
  \item{wrapper.el}{ Wrapper elements; generally set by \code{\link{MARSS}} }
  \item{allowed}{ Allowed constraints.  This changes depending on the fitting method the user has specified (in \code{\link{MARSS}} or in the \code{\link{marssMLE}} object).  Lists of for \code{allowed} for different fitting methods is set in \code{MARSSsettings.R}. Use \code{MARSS:::allowed} to view since \code{allowed} is hidden.}
  \item{silent}{ Suppresses errors and warnings printing. }
}
\details{
  Called by \code{\link{popWrap}} to ensure that user inputs are valid and can be handled by \code{\link{as.marssm}}. 
}
\value{
  TRUE if object passes all checks.
}
\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
} 
\seealso{ \code{\link{popWrap}} \code{\link{as.marssm}} }
\examples{
  \dontrun{
  ## Error:
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  wrapperObj = popWrap(dat, allowed=MARSS:::allowed$kem, constraint=list(Z="wrong"))
  }
}

