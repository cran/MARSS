\name{MARSSoptim}
\alias{MARSSoptim}
\alias{neglogLik}
\title{Parameter estimation for MARSS models using optim}
\description{
  Parameter estimation for MARSS models using R's \code{\link{optim}} function.   This allows access to R's quasi-Newton algorithms available via the \code{\link{optim}} function.  The \code{MARSSoptim} is called when \code{\link{MARSS}} is called with \code{method="BFGS"}.  Only diagonal Q and R matrices are allowed, if they are estimated.  This is a base function in the \code{\link{MARSS-package}}.  \code{neglogLik} is a helper function for \code{MARSSoptim} that returns the negative log-likelihood given a vector of the estimated parameters and a \code{\link{marssMLE}} object.
}
\usage{
MARSSoptim(MLEobj)
neglogLik(x, MLEobj)
}
\arguments{
  \item{MLEobj}{An object of class \code{\link{marssMLE}}.}
  \item{x}{An vector of the estimated parameters as output by \code{\link{MARSSvectorizeparam}}.}
}
\details{
Objects of class \code{\link{marssMLE}} may be built from scratch but are easier to construct using \code{\link{MARSS}} with \code{MARSS(..., fit=FALSE, method="BFGS")}.

  Options for \code{\link{optim}} are passed in using \code{MLEobj$control}.  See \code{\link{optim}} for a list of that function's control options.  If \code{lower} and \code{upper} for \code{\link{optim}} need to be passed in, they should be passed in as part of \code{control} as \code{control$lower} and \code{control$upper}. Additional \code{control} arguments affect printing and initial conditions.
  \describe{
      \item{\code{MLEobj$control$iter.V0}}{ This is the value to which V0 will be set during the maximization algorithm (must not be equal to 0). See Details.}
      \item{\code{MLEobj$control$MCInit}}{ If TRUE, Monte Carlo initialization will be performed by \code{\link{MARSSmcinit}}. } 
      \item{\code{MLEobj$control$numInits}}{ Number of random initial value draws to be used with \code{\link{MARSSmcinit}}. Ignored if \code{control$MCInit=FALSE}. } 
      \item{\code{MLEobj$control$numInitSteps}}{ Maximum number of EM iterations for each random initial value draw to be used with \code{\link{MARSSmcinit}}. Ignored if \code{control$MCInit=FALSE}. }
      \item{\code{MLEobj$control$boundsInits}}{ Length 6 list. Each component is a length 2 vector of bounds on the uniform distributions from which initial values will be drawn to be used with \code{MARSSmcinit()}. Ignored if \code{control$MCInit=FALSE}. See Examples.}
      \item{\code{MLEobj$control$silent}}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }    
    }
}
\value{
  The \code{\link{marssMLE}} object which was passed in, with additional components:
  \item{method}{String "BFGS".}
  \item{kf}{Kalman filter output. }
  \item{iter.record}{If \code{MLEobj$control$trace = TRUE}, then this is the \code{$message} value from \code{\link{optim}}.}
  \item{numIter}{Number of iterations needed for convergence.}
  \item{convergence}{ Did estimation converge successfully? 
    \describe{
      \item{convergence=0}{ Converged in less than \code{MLEobj$control$maxit} iterations and no evidence of degenerate solution. } 
      \item{convergence=1}{ Maximum number of iterations \code{MLEobj$control$maxit} was reached before \code{MLEobj$control$abstol} condition was satisfied. } 
      \item{convergence=10}{ Some of the variance elements appear to be degenerate.  T} 
      \item{convergence=52}{ The algorithm was abandoned due to errors from the "L-BFGS-B" method.  } 
      \item{convergence=53}{ The algorithm was abandoned due to numerical errors in the likelihood calculation from \code{\link{MARSSkf}}.  This happens frequently with "BFGS" and can sometimes be helped with a better initial condition.  Try using the Kalman-EM algorithm first (\code{method="kem"}), and then using the parameter estimates from that to as initial conditions for \code{method="BFGS"}.  } 
      }
  }
  \item{logLik}{Log-likelihood.}
  \item{states}{State estimates from the Kalman filter.}
  \item{states.se}{Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) Shumway & Stoffer.}
  \item{errors}{Any error messages.}
}
\section{Discussion}{
 The function only returns parameter estimates.  To compute CIs, use \code{\link{MARSSparamCIs}} but if you use parametric or non-parametric bootstrapping with this function, it will use the Kalman-EM algorithm to compute the bootstrap parameter estimates!  The quasi-Newton estimates are too fragile for the bootstrap routine since one often needs to search to find a set of initial conditions that work (i.e. don't lead to numerical errors).
 
 Estimates from \code{MARSSoptim} (which come from \code{\link{optim}}) should be checked against estimates from Kalman-EM algorithm.  If the quasi-Newton algorithm works, it will tend to find parameters with higher likelihood faster than the Kalman-EM algorithm.  However, the MARSS likelihood surface can be multimodal with sharp peaks at degenerate solutions where a Q or R diagonal element equals 0.  The quasi-Newton algorithm tends to find and gets stuck on these peaks even when they are not the maximum.  Neither an initial conditions search nor starting near the known maximum (or from the parameters estimates after the Kalman-EM algorithm) will necessarily solve this problem.   Thus it is wise to check against Kalman-EM estimates to ensure that the BFGS estimates are close to the MLE estimates.  
 
 Note this is mainly a problem if the time series are short or very gappy.  If the time series are long, then the likelihood surface should be nice with a single interior peak.  In this case, the quasi-Newton algorithm works well but it can still be sensitive (and slow) if not started with a good initial condition.  Thus starting it with the estimates from the Kalman-EM algorithm is often desirable.
 
 The V0 matrix is set to a diffuse prior if x0 is estimated (in which case V0 must be zero to treat x0 as fixed but unknown).  V0 is reset to zero when the final likelihood and state estimates are computed via \code{\link{MARSSkf}}.  See discussion in the manual and \code{\link{MARSSkem}}. \code{MARSSoptim} only allows diagonal Q and R matrices, if they are estimated. 
}
\author{
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{  \code{\link{MARSS}} \code{\link{MARSSkem}}  \code{\link{marssMLE}}  \code{\link{optim}}
}
\examples{
dat = t(harborSealWA)
dat = dat[2:nrow(dat),] #remove the year row
#fit a model with 1 hidden state where obs errors are iid
bfgsfit = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1))), 
   method="BFGS")

#fit a model with Kalman-EM and then use that fit as the start for BFGS
kemfit = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1))))
bfgsfit = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1))), 
   inits=kemfit$par, method="BFGS")
}
