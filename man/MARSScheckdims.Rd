\name{MARSScheckdims}
\alias{MARSScheckdims}
\alias{MARSScheckpar}
\title{MARSS Utilities}
\description{
  This is a utility function in the \code{\link{MARSS-package}} for checking MARSS matrix dimensions and parameter lists.
}
\usage{
MARSScheckdims(el, target, n, m)
MARSScheckpar(parList, n, m)
}
\arguments{
  \item{el}{ Name of a MARSS list element ("A", "B", "Q", "R", "U", "V0", "x0", "Z"). }
  \item{target}{ List to be checked. }
  \item{n}{ Number of time series. }
  \item{m}{ Number of state processes. }
  \item{parList}{ MARSS parameter list. }
}
\value{
  TRUE if no problems; otherwise a message describing the problems.
}
\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
}


