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
Passing in \code{form="dfa"} is invokes a helper function to create that model and creates the Z matrix for the user.  Q is by definition identity, x0 is zero and V0 is diagonal with large variance (5). U is zero, A is zero, and covariates only enter the Y equation.  Because U and A are 0, the data should have mean 0 (de-meaned) otherwise one is likely to be creating a structurally inadequate model (i.e. the model implies that the data have mean = 0, yet data do not have mean = 0 ).
}
\usage{
MARSS.dfa(MARSS.inputs)
}

\arguments{
  \item{MARSS.inputs}{ A list of the inputs to a MARSS() call. When \code{form="dfa"} in addition to the usual \code{\link{MARSS}} inputs, the user can also pass in:
      \itemize{
    \item \code{demean}  Logical.  Default is TRUE, which means the data will be demeaned.
    \item \code{z.score}  Logical.  Default is TRUE, which means the data will be z-scored (demeaned and variance standardized to 1).
    }
    The \code{model} argument of the \code{\link{MARSS}} call is constrained so that the DFA model above is specified and allows an additional element, \code{m}, that allows the user to specify the number of hidden state variables.  It is not necessariy for the user to specify \code{Z} as the helper function will create a \code{Z} appropriate for a DFA model.
  }
}
\details{
The user can pass in altername forms for some of the model parameters.  In addition, m can be passed in to set the number of state variables. The following are allowed (and not allowed):
    \itemize{
    \item \code{B}  Cannot be changed or passed in via model argument.
    \item \code{U}  Cannot be changed or passed in via model argument.
    \item \code{Q}  Cannot be changed or passed in via model argument.
    \item \code{Z}  Can be passed in a matrix if the user does not want a default DFA \code{Z} matrix.  There are many equivalent ways to construct a DFA \code{Z} matrix.  We use Zuur et al.'s form.
    \item \code{A}  Can be "unequal", "zero" or a matrix.
    \item \code{R}  Can be "identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or a matrix
    \item \code{x0} Can be "unconstrained", "unequal", "zero", or a matrix.
    \item \code{V0} Can be "identity", "zero", or a matrix.
    \item \code{tinitx} Can be 0 or 1.
    \item \code{m} Can be 1 to n (the number of y time-series).
    }

See the User Guide chapter on Dynamic Factor Analysis for examples of of using \code{form="dfa"}.
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
