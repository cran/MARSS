\name{MARSS-package}
\alias{MARSS-package}
\docType{package}
\title{
Multivariate Autoregressive State-Space Model Estimation
}
\description{
The MARSS package fits constrained and unconstrained multivariate autoregressive time-series models to multivariate time series data.  To open the user guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}.  To open an overview page with package information, type \code{RShowDoc("index",package="MARSS")}.    

The MARSS model is
\describe{
\item{x(t+1) = B x(t) + U + w(t), where w(t) ~ MVN(0,Q)}{}
\item{y(t) = Z x(t) + A + v(t), where v(t) ~ MVN(0,R)}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0) }{} 
}
The parameters, hidden state processes (x), and observations (y) are matrices:
  \itemize{
    \item x(t) is m x 1 
    \item y(t) is n x 1 (m<=n)
    \item Z is n x m 
    \item B is m x m 
    \item U is m x 1 
    \item Q is m x m 
    \item A is n x 1 
    \item R is n x n 
    \item x0 is m x 1 
    \item V0 is m x m 
  }
The package functions estimate the model parameters using an EM algorithm (primarily but see \code{\link{MARSSoptim}}).  Parameters may be constrained to have shared elements (elements which are constrained to have the same value) or fixed elements (with the other elements estimated).  The states and smoothed state estimates are provided via a Kalman filter and smoother.  Bootstrapping, confidence interval estimation, bias estimation, model selection and simulation functions are provided.   The main user interface to the package is the top-level function \code{\link{MARSS}}.  
}
\details{
Important MARSS functions:
\describe{
      \item{\code{\link{MARSS}}}{ Top-level function for specifying and fitting MARSS models. } 
      \item{\code{\link{MARSSsimulate}}}{ Produces simulated data from a MARSS model. } 
      \item{\code{\link{MARSSkem}}}{ Estimates MARSS parameters using an EM algorithm. } 
      \item{\code{\link{MARSSkf}} and \code{\link{MARSSkfas}}}{ Kalman filter and smoother. } 
      \item{\code{\link{MARSSoptim}}}{ Estimates MARSS parameters using a quasi-Newton algorithm via \code{\link{optim}}. } 
      \item{\code{\link{MARSSaic}}}{ Calculates AICc, AICc, and various bootstrap AICs. } 
      \item{\code{\link{MARSSboot}}}{ Creates bootstrap MARSS parameter estimates. }
      \item{\code{\link{MARSSparamCIs}}}{ Computes confidence intervals for maximum-likelihood estimates of MARSS parameters. }    
    }
}
\author{
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov, 
  
  kellie(dot)wills(at)noaa(dot)gov
}
\references{
The MARSS user guide:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112.         

Type \code{RShowDoc("UserGuide",package="MARSS")} at the R command line to open the MARSS user guide.

Type \code{RShowDoc("index",package="MARSS")} to see all the package documentation, tutorials, and case study scripts.
}

\keyword{ package }
