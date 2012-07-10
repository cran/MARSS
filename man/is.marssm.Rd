\name{is.marssm}
\alias{is.marssm}
\alias{as.marssm}
\alias{marssm}
\title{ Model Objects }
\description{
These are model objects and utility functions for model objects in the package \code{\link{MARSS-package}}.  
Users would not normally work directly with these functions.
  \code{is.marssm()} ensures model consistency.
  \code{MARSS.formname()} translates a model list as passed in call to MARSS() into a marssm model object.
}
\usage{
is.marssm(modelObj)
}

\arguments{
  \item{modelObj}{ An object of class marssm. }
}
\details{
A \code{marssm} object is an R representation of a MARSS model along with the data.  
    Data for \code{marssm()} consists of multivariate time series data in which time is across columns and the n observed time series are in the n different rows.

The MARSS model is
\describe{
\item{x(t+1) = B(t) x(t) + U(t) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + A(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(x0, V0)}{} 
}

The marssm object describes this MARSS model but written in vec form:
\describe{
\item{x(t+1) = kron(x(t),I)(f_b(t)+D_b(t)b) + (f_u(t)+D_u(t)u) + w(t), where w(t) ~ MVN(0,Q)}{}
\item{vec(Q) = f_q(t)+D_q(t)q}{} 
\item{y(t) = kron(x(t),I)(f_z(t)+D_z(t)z) + (f_a(t)+D_a(t)a) + v(t), where v(t) ~ MVN(0,R)}{}
\item{vec(R) = f_r(t)+D_r(t)r}{} 
\item{x(1) ~ MVN(f_p+D_p p, V0)}{}
\item{vec(V0) = f_l+D_l l}{}  
}
In the marssm object, f(t) + D(t)m, is the vec of a matrix M(t), so f_b(t)+D_b(t)b would be vec(B(t)). The estimated parameters are in the column vectors: b, u, q, z, a, r, p, and l.  Each matrix M(t) is f(t)+D(t)m so is the sum of a fixed part f(t) and the linear combination, D(t), of the free (or estimated) part m.  

The vec form of the MARSS model is specified by 3D matrices for each f and D for each parameter: B, U, Q, Z, A, R, x0, V0.  The number of columns in the D matrix for a parameter determines the number of estimated values for that parameter.  
The first dimension for f (\code{fixed}) and D (\code{free}) must be:
  \describe{
    \item{Z}{ n x m }
    \item{B, Q, and V0}{ m x m }
    \item{U and x0}{ m x 1 }
    \item{A}{ n x 1 }
    \item{R}{ n x n }
  }
  The third dimension of f (\code{fixed}) and D (\code{free}) is either 1 (if not time-varying) or TT (if time-varying).   The second dimension of f (\code{fixed}) is always 1, while the second dimension of D  (\code{free}) depends on how many values are being estimated for a matrix.  It can be 0 (if the matrix is fixed) or up to the size of the matrix (if all elements are being estimated).

Many different types of multivariate time series models can be rewritten in this form.  Users however are not expected to write their MARSS model in this form.  Instead, the main \code{\link{MARSS}} function calls helper functions of the form \code{MARSS.formname} to transform user inputs written in forms more familiar to the user into the vec form needed for the algorithms in the MARSS package.
See the user guide (\code{RShowDoc("UserGuide",package="MARSS")}) for many examples.
}
\value{
  An object of class "marssm".
  \item{data}{ Data supplied by user. }
  \item{fixed}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{free}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{X.names}{ A m length vector of names for the X states.} 
  \item{miss.value}{ Specifies missing value representation. Default is NA}
  \item{tinitx}{ Specifies whether the initial state (specified with x0 and V0) is at t=0 or t=1. }
}
\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
}

\seealso{ \code{\link{MARSS}}, \code{\link{MARSS.marxss}} } 


