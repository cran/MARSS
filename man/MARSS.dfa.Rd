\name{MARSS.dfa}
\alias{MARSS.dfa}
\title{ Multivariate Dynamic Factor Analysis }
\description{
The Dynamic Factor Analysis model in MARSS is
\describe{
\item{x(t+1) = x(t) + w(t), where w(t) ~ MVN(0,I)}{}
\item{y(t) = Z(t) x(t) + D(t)d(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(0, 5*I)}{} 
}
}
\usage{
MARSS.dfa(MARSS.inputs)
}

\arguments{
  \item{MARSS.inputs}{ A list of the inputs to a MARSS() call.}
}
\details{
See the User Guide chapter on Dynamic Factor Analysis for a description of the MARSS.dfa function.
}
\value{
  An object of class "marssm".
  \item{data}{ Data supplied by user. }
  \item{fixed}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{free}{ A list with 8 matrices Z, A, R, B, U, Q, x0, V0. }
  \item{X.names}{ A m length vector of names for the X states.} 
  \item{miss.value}{ Specifies missing value representation. Default is NA}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA. 
  }

\seealso{ \code{\link{MARSS}} } 
