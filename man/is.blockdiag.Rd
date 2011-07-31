\name{is.blockdiag}
\alias{is.blockdiag}
\alias{is.blockequaltri}
\alias{is.blockunconst}
\alias{is.diagonal}
\alias{is.equaltri}
\alias{makediag}
\alias{takediag}
\alias{is.design}
\alias{is.fixed}
\alias{vec}
\alias{unvec}
\alias{is.wholenumber}
\alias{as.design}
\alias{Imat}
\alias{is.identity}
\title{ Matrix Utilities }
\description{
  Matrix utilities for MARSS functions in the \code{\link{MARSS-package}}.
}
\usage{
is.blockdiag(x)
is.blockequaltri(x, uniqueblocks=FALSE)
is.blockunconst(x, uniqueblocks=FALSE)
is.diagonal(x, na.rm=FALSE)
is.equaltri(x)
makediag(x, nrow=NA)
takediag(x)
is.design(x)
is.fixed(x)
is.identity(x)
vec(x)
unvec(x, dim=NULL)
is.wholenumber(x, tol = .Machine$double.eps^0.5)
as.design(fixed,free)
Imat(x)
}
\arguments{
  \item{x}{ A matrix (or vector for '\code{makediag}'). }
  \item{na.rm}{ How to treat NAs in the block diag test. }
  \item{dim}{ Matrix dimensions. }
  \item{fixed}{ A fixed matrix per the MARSS specification for fixed matrix syntax. }
  \item{free}{ A free matrix per the MARSS specification for free matrix syntax. }
  \item{nrow}{ Number of rows. }
  \item{tol}{ Tolerance. }
  \item{uniqueblocks}{ Must blocks be unique? }
}
\details{
  '\code{is...}' tests for various matrix properties.
  \code{vec(x)} creates a column vector from a matrix per the standard vec math function.
  \code{unvec(c,dim)} takes the vector c and creates a matrix with the specified dimensions.
  \code{as.design(fixed,free)} returns the fixed vector and design matrix for a fixed/free pair.
  \code{Imat(nrow)} returns the identity matrix of dimension nrow.
}
\value{
  '\code{makediag(x)}': a matrix with diagonal x.
  '\code{takediag(x)}': the diagonal from matrix x.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}

