\name{MARSSkem}
\alias{MARSSkem}
\title{ Maximum Likelihood Estimation for Multivariate Autoregressive State-Space Models }
\description{
  \code{MARSSkem()} performs maximum-likelihood estimation, using an EM algorithm for constrained and unconstrained MARSS models.  This is one of the base functions in the \code{\link{MARSS-package}}.
}
\usage{
MARSSkem(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
}
\details{
Objects of class \code{\link{marssMLE}} may be built from scatch but are easier to construct using \code{\link{MARSS}} with \code{MARSS(..., fit=FALSE)}.

  Options for \code{MARSSkem()} may be set using \code{MLEobj$control}, as follows:
  \describe{
      \item{\code{MLEobj$control$minit}}{ Minimum number of EM iterations. You can use this to force the algorithm to do a certain number of iterations.  This is helpful if your soln is not converging.} 
      \item{\code{MLEobj$control$maxit}}{ Maximum number of EM iterations. } 
      \item{\code{MLEobj$control$abstol}}{ Tolerance for log-likelihood change.  If log-likelihood changes less than this amount relative to the previous iteration, the EM algorithm exits. } 
      \item{\code{MLEobj$control$iter.V0}}{ This is the value to which V0 will be set during the EM algorithm (!=0). See Details.}
      \item{\code{MLEobj$control$trace}}{ A positive integer. If not 0, a record will be created of each variable over all EM iterations and detailed warning messages (if appropriate) will be printed. }
      \item{\code{MLEobj$control$safe}}{ If TRUE, \code{MARSSkem} will rerun \code{\link{MARSSkf}} after each individual parameter update rather than only after all parameters are updated.  The latter is slower and unnecessary for many models, but in some cases, the safer and slower algorithm is needed because the ML parameter matrices have high condition numbers. }
      \item{\code{MLEobj$control$MCInit}}{ If TRUE, Monte Carlo initialization will be performed by \code{\link{MARSSmcinit}}. } 
      \item{\code{MLEobj$control$numInits}}{ Number of random initial value draws to be used with \code{\link{MARSSmcinit}}. Ignored if \code{control$MCInit=FALSE}. } 
      \item{\code{MLEobj$control$numInitSteps}}{ Maximum number of EM iterations for each random initial value draw to be used with \code{\link{MARSSmcinit}}. Ignored if \code{control$MCInit=FALSE}. }
      \item{\code{MLEobj$control$boundsInits}}{ Length 6 list. Each component is a length 2 vector of bounds on the uniform distributions from which initial values will be drawn to be used with \code{MARSSmcinit()}. Ignored if \code{control$MCInit=FALSE}. See Examples.}
      \item{\code{MLEobj$control$silent}}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }    
    }
}
\value{
  The \code{\link{marssMLE}} object which was passed in, with additional components:
  \item{method}{String "kem".}
  \item{kf}{Kalman filter output. }
  \item{iter.record}{If \code{MLEobj$control$trace = TRUE}, a list with \code{par} = a record of each estimated parameter over all EM iterations and \code{logLik} = a record of the log likelikelihood with V0 set to the value of \code{iter.V0}.  Note this is different than the log likelihood with V0 = 0 (which is the final log likelihood value).}
  \item{numIter}{Number of iterations needed for convergence.}
  \item{convergence}{ Did estimation converge successfully? 
    \describe{
      \item{convergence=0}{Converged in less than \code{MLEobj$control$maxit} iterations and no evidence of degenerate solution.}  
      \item{convergence=1}{Maximum number of iterations \code{MLEobj$control$maxit} was reached before \code{MLEobj$control$abstol} condition was satisfied.} 
      \item{convergence=10}{Some of the variance elements appear to be degenerate.  The \code{MLEobj$control$abstol} condition was satisfied 
          before the maximum number of iterations \code{MLEobj$control$maxit} was reached, but some of the variance elements appear to still be changing. } 
      \item{convergence=52}{The algorithm was abandoned due to numerical errors.  Usually this means one of the variances either went to zero or to all elements being equal.  Try setting \code{control$trace=1} to view a detailed error report.}  
      }
  }
  \item{logLik}{Log-likelihood.}
  \item{states}{State estimates from the Kalman filter.}
  \item{states.se}{Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) Shumway & Stoffer.}
  \item{errors}{Any error messages.}
}

\section{Discussion}{
 To ensure that the global maximum-likelihood values are found, it is recommended that initial parameter values be set using Monte Carlo initialization (\code{MLEobj$control$MCInit} = TRUE), particularly if the model is not a good fit to the data.  This requires more compuation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.   For many models, this is unnecessary, but answers should be checked using an initial conditions search before reporting final values.

 \code{MARSSkem()} calls a Kalman filter/smoother (\code{\link{MARSSkf}}) for hidden state estimation.   The algorithm allows two options for the initial state conditions: x at t=1.  Either x0 is treated as fixed but unknown (estimated); in this case, \code{fixed$V0=0} and x0 is estimated.  This is the default behavior.  In the second case, the initial conditions are specified with a known prior, \code{fixed$x0} and \code{fixed$V0}, and x0 is not estimated.   In the first case, x0 fixed but unknown and to be estimated, V0=0 but the algorithm cannot be run with this since x0 would never move from its initial value at iteration=1.  Instead, during the EM iterations, a diffuse prior is used.  This is done by setting \code{MLEobj$control$iter.V0} to a large value, say 10.  A small value will cause the algorithm to converge very slowly and 0 will generate an error. Before reporting the final log likelihood, \code{MARSSkem()} runs the Kalman filter with the maximum-likelihood parameter estimates (using the diffuse prior) and V0 = 0 to obtain the correct likelihood.  This approach works well and is easy to implement.  However when one fits a model with x0 having shared values AND an unconstrained R matrix, an ill-conditioned matrix tends to appear in one of the steps of the Kalman filter algorithm (because V0 must have elements with 100 percent correlation if you say that some x0's have the same value).  \code{MARSSkem()} will report warnings and errors if this happens.  Switching to a diagonal R or an unconstrained x0 will fix the ill-conditioning problems.  See the manual for discussion about how V0 is treated in the algorithm.
 
 If you get errors, it generally means that the solution involves an ill-conditioned matrix.    For example, your Q or R matrix is going to a value in which all elements have the same value, for example zero.  If for example, you tried to fit a model with fixed and high R matrix and the variance in that R matrix was much higher than what is actually in the data, then you might drive Q to zero.   Also if you try to fit a structurally inadequate model, then it is not unusual that Q will be driven to zero.  For example, if you fit a model with 1 hidden state trajectory to data that clearly have 2 quite different hidden state trajectories, you might have this problem.  Comparing the likelihood of this model to a model with more structural flexibility should reveal that the structually inflexible model is inadequate (much lower likelihood).
 
 If you get a warning that the solution appears to be degenerate, it means that some of the elements in Q or R are going to zero and the log-likelihood is changing very slowly. 
 You can either try decreasing \code{control$abstol} dramatically (e.g. 1E-6), use a Newton finisher \code{\link{MARSSoptim}}\code{(MLEobj)}, or fix
 the degenerate values to something very small (e.g. 1E-8) and re-estimate.  Try \code{\link{find.degenerate}}(MLEobj) using the output from the \code{MARSSkem}
 call to find the degenerate elements.
     }

\references{ 
  R. H. Shumway and D. S. Stoffer (2006).  Chapter 6 in Time Series Analysis and its Applications.  Springer-Verlag, New York.
  
  Ghahramani, Z. and Hinton, G. E. (1996) Parameter estimation for linear dynamical systems. Technical Report CRG-TR-96-2, University of
Totronto, Dept. of Computer Science. 

Harvey, A. C. (1989) Chapter 5 in Forecasting, structural time series models and the
Kalman filter. Cambridge University Press, Cambridge, UK.
  
  The user manual:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   type \code{show.doc(MARSS, manual)} to see.          

  The EM algorithm: Holmes, E. E. (2010).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models.  type \code{show.doc(MARSS, Holmes2010)} to see.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{MARSSmcinit}}, \code{\link{MARSSkf}}, \code{\link{marssMLE}}, \code{\link{MARSSoptim}},  \code{\link{find.degenerate}}
  }  
\examples{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),] 
  #you can use MARSS to construct a proper marssMLE object.
  MLEobj = MARSS(dat, fit=FALSE) 
  kemfit = MARSSkem(MLEobj)
  \dontrun{
  #see what a Newton method would find; 
  # wrapped in try because it tends to throw numerical errors
  bfgsfit=kemfit; bfgsfit$method="BFGS"
  bfgsfit = try(MARSSoptim(bfgsfit))
  #look for degenerate estimates
  #Use MARSS to do the fit
  kemfit2 = MARSS(dat, control=list(minit=100)) 
  #this will make a plot that will show if the vars converged.
  find.degenerate(kemfit2) 
  }
}


