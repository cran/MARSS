\name{ MARSSoptions }
\alias{MARSSoptions}
\title{ Change MARSS Defaults Utility }
\description{
  Utility in the \code{\link{MARSS-package}} to change the defaults, including default model structure, for the MARSS function.
}
\usage{
MARSSoptions(..., method="kem")
}
\arguments{
  \item{...}{ A name or list of names in \code{alldefaults}.  To see what  \code{alldefaults} looks like, type \code{MARSSoptions()}. This is the same as \code{(MARSS:::alldefaults)[[method]]}}.
  \item{method}{ Estimation method. MARSS 1.0 allows method="kem" and "BFGS". }
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
\dontrun{
#Change the defaults maxit value
MARSSoptions(control=list(maxit=5000))

#Change the defaults minit, maxit and abstol value; setting abstol means the program defaults to using the delta logLik convergence test
MARSSoptions(control=list(minit=100, maxit=5000, abstol=0.01))

#Change lots of different defaults
MARSSoptions(control=list(minit=100, maxit=5000, abstol=0.01), constraint=list(U="unequal",Q="unconstrained",R=diag(0.01,4)))

#Show the control defaults
MARSSoptions("control")

#Show all the defaults for method="BFGS"; this doesn't change the default method
MARSSoptions(method="BFGS")
}
}

