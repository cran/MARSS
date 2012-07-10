\name{is.blockdiag}
\alias{makediag}
\alias{takediag}
\alias{is.blockdiag}
\alias{is.identity}
\alias{is.diagonal}
\alias{is.equaltri}
\alias{is.blockequaltri}
\alias{is.blockunconst}
\alias{is.design}
\alias{is.fixed}
\alias{vec}
\alias{unvec}
\alias{is.wholenumber}
\alias{Imat}
\alias{rwishart}
\alias{is.identity}
\alias{ mystrsplit }
\alias{ convert.model.mat }
\alias{fixed.free.to.formula }
\alias{ matrix.power }
\alias{ sub3D }
\alias{ pinv }
\alias{ pcholinv }
\alias{ is.solvable }

\title{ Matrix Utilities }
\description{
  Matrix utilities for MARSS functions in the \code{\link{MARSS-package}}.  These are not exported but can be accessed using the \code{MARSS:::} prefix.
}
\usage{
is.blockdiag(x)
is.blockequaltri(x, uniqueblocks=FALSE)
is.blockunconst(x, uniqueblocks=FALSE)
is.identity(x, dim=NULL)
is.diagonal(x, na.rm=FALSE)
is.equaltri(x)
makediag(x, nrow=NA)
takediag(x)
is.design(x, strict=TRUE, dim=NULL, zero.rows.ok=FALSE, zero.cols.ok=FALSE)
is.fixed(x, by.row=FALSE)
is.identity(x, dim=NULL)
vec(x)
unvec(x, dim=NULL)
is.wholenumber(x, tol = .Machine$double.eps^0.5)
Imat(x)
rwishart(nu, V)
mystrsplit(x)
convert.model.mat(param.matrix)
fixed.free.to.formula(fixed,free,dim)
matrix.power(x, n)
sub3D(x,dim1,dim2,t=1)
pinv(x)
pcholinv(x)
is.solvable(A,y=NULL)
}
\arguments{
  \item{x, A, y}{ A matrix (or vector for '\code{makediag}' or string for '\code{mystrsplit}'.). }
  \item{na.rm}{ How to treat NAs in the block diag test. }
  \item{dim}{ Matrix dimensions. Some functions will take the vec of a matrix.  In this case, the optional dim arg specifies the matrix dimensions. }
  \item{fixed}{ A fixed matrix per the MARSS specification for fixed matrix syntax. }
  \item{free}{ A free matrix per the MARSS specification for free matrix syntax. }
  \item{nrow}{ Number of rows. }
  \item{tol}{ Tolerance. }
  \item{uniqueblocks}{ Must blocks be unique? }
  \item{t}{ The time index or third dimension of a 3D matrix }
  \item{nu, V}{ Parameters of a Wishart distribution. }
  \item{param.matrix}{ The list matrix version of a time-invariant MARSS model. }
  \item{dim1, dim2}{ The first and second dimensions of a 3D matrix. }
  \item{n}{ An interger for the power function. }
  \item{ zero.rows.ok, zero.cols.ok }{ Means the design matrix can have all zero rows or columns. }
  \item{ strict }{ Specifies whether the design matrix must be only 0s and 1s. }
  \item{ by.row }{ For is.fixed, reports whether is.fixed by row rather than for the whole matrix. }
}
\details{
  '\code{is...}' tests for various matrix properties.
  \code{vec(x)} creates a column vector from a matrix per the standard vec math function.
  \code{unvec(c,dim)} takes the vector c and creates a matrix with the specified dimensions.
  \code{Imat(nrow)} returns the identity matrix of dimension nrow.
  \code{ fixed.free.to.formula } takes a fixed and free pair and constructs a list matrix with formulas in each matrix element.
  \code{ convert.model.mat } takes a list matrix with formulas in each element and converts to a fixed/free pair.
  \code{ sub3D } returns a 2D matrix after subsetting a 3D matrix on the third (time) dimension.
  \code{ mystrsplit } is a customized string splitter used by \code{convert.model.mat}.
  \code{ rwishart } generates random draws from a wishart distribution.
  \code{ matrix.power } is a faster way to get the n-th power of a matrix.
  \code{ pinv } is the pseudoinverse based on singular value decomposition PInv=UD^+V' where a diagonal matrix with non-zero diagonal values of D (from svd) replaced with 1/D.
  \code{ pcholinv } is the inverse based on the Cholsky decomposition but modified to allow 0s on the diagonal of x (with corresponding 0 row/column).  These appear as 0 row/columns in the returned inverse.
  \code{ is.solvable } returns information on the solvability of the linear system y=Ax using the SVD decomposition.
}
\value{
  '\code{makediag(x)}': a matrix with diagonal x.
  '\code{takediag(x)}': the diagonal from matrix x.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}

