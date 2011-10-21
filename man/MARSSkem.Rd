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

  Options for \code{MARSSkem()} may be set using \code{MLEobj$control}.  The commonly used elements of \code{control} are follows (see \code{\link{marssMLE}}:
  \describe{
      \item{\code{minit}}{ Minimum number of EM iterations. You can use this to force the algorithm to do a certain number of iterations.  This is helpful if your soln is not converging.} 
      \item{\code{maxit}}{ Maximum number of EM iterations. } 
      \item{\code{min.iter.conv.test}}{ The minimum number of iterations before the log-log convergence test will be computed.  If \code{maxit} is set less than this, then convergence will not be computed (and the algorithm will just run for maxit iterations). } 
      \item{\code{kf.x0}}{ Whether to set the prior at t=0 ("x00") or at t=1 ("x10"). The default is "x00". } 
      \item{\code{conv.test.deltaT}}{ The number of iterations to use in the log-log convergence test.  This defaults to 9. } 
      \item{\code{abstol}}{ Tolerance for log-likelihood change for the delta logLik convergence test.  If log-likelihood changes less than this amount relative to the previous iteration, the EM algorithm exits. This is normally (default) set to NULL and the log-log convergence test is used instead.} 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{ A positive integer. If not 0, a record will be created of each variable over all EM iterations and detailed warning messages (if appropriate) will be printed. }
      \item{\code{safe}}{ If TRUE, \code{MARSSkem} will rerun \code{\link{MARSSkf}} after each individual parameter update rather than only after all parameters are updated.  The latter is slower and unnecessary for many models, but in some cases, the safer and slower algorithm is needed because the ML parameter matrices have high condition numbers. }
      \item{\code{MCInit}}{ If TRUE, Monte Carlo initialization will be performed by \code{\link{MARSSmcinit}}. } 
      \item{\code{numInits}}{ Number of random initial value draws to be used with \code{\link{MARSSmcinit}}. Ignored if \code{MCInit=FALSE}. } 
      \item{\code{numInitSteps}}{ Maximum number of EM iterations for each random initial value draw to be used with \code{\link{MARSSmcinit}}. Ignored if \code{MCInit=FALSE}. }
      \item{\code{boundsInits}}{ Length 6 list. Each component is a length 2 vector of bounds on the uniform distributions (for B, U, and A) from which initial values will be drawn to be used with \code{MARSSmcinit()}. For R and Q, the 2 bounds specify the df and S(=diag(bound[2],m)) for a wishart distribution. Ignored if \code{MCInit=FALSE}. See Examples.}
      \item{\code{silent}}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }    
    }
}
\value{
  The \code{\link{marssMLE}} object which was passed in, with additional components:
  \item{method}{String "kem".}
  \item{kf}{Kalman filter output. }
  \item{iter.record}{If \code{MLEobj$control$trace = TRUE}, a list with \code{par} = a record of each estimated parameter over all EM iterations and \code{logLik} = a record of the log likelikelihood at each iteration. }
  \item{numIter}{Number of iterations needed for convergence.}
  \item{convergence}{ Did estimation converge successfully? 
    \describe{
      \item{convergence=0}{Converged in less than \code{MLEobj$control$maxit} iterations and no evidence of degenerate solution using the log-log plot test.}  
      \item{convergence=1}{Maximum number of iterations \code{MLEobj$control$maxit} was reached before \code{MLEobj$control$abstol} condition was satisfied. This value can only be output if \code{MLEobj$control$abstol} is passed in.  If the default log-log convergence test is being used, convergence=10 is returned when maxit is reached.}
      \item{convergence=2}{No convergence diagnostic were computed because there are problems with the MLE object.}  
      \item{convergence=3}{No convergence diagnostic were computed because \code{MLEobj$control$maxit} was set to less than \code{control$min.iter.conv.test}.  This is not an error.}  
      \item{convergence=10}{Some of the parameter estimates did not converge (based on the log-log plot test) before \code{MLEobj$control$maxit} was reached. This is not an error per se. Degenerate solutions are fine but the MARSS algorithm will not compute the proper likelihood for a degenerate solution. }
      \item{convergence=11}{Some of the parameter estimates did not converge (based on the log-log plot test) even though \code{MLEobj$control$abstol} was reached. This value can only be output if \code{MLEobj$control$abstol} is passed in. Try not setting \code{MLEobj$control$abstol} so that the log-log convergence test is used instead.  If it takes too long, try making  \code{MLEobj$control$maxit} smaller (like 500). }
      \item{convergence=12}{\code{MLEobj$control$abstol} was reached but the log-log plot test returned NAs.  This is an odd error and you should set \code{control$trace=TRUE} and look at the outputted \code{$iter.record} to see what is wrong. This value can only be output if \code{MLEobj$control$abstol} is passed in. }
      \item{convergence=52}{The algorithm was abandoned due to numerical errors.  Usually this means one of the variances either went to zero or to all elements being equal.  This is not an error per se.  Most likely it means that your model is not very good for your data (too inflexible or too many parameters).  Try setting \code{control$trace=1} to view a detailed error report.}
      \item{convergence=62}{The algorithm was abandoned due to errors in the log-log convergence test.  You should not get this error (it is included for debugging purposes to catch improper arguments passed into the log-log convergence test).  If you get it, try passing in \code{control$abstol} to use delta logLik as the convergence test.}
      \item{convergence=72}{Other convergence errors.  This is included for debugging purposes to catch misc. errors.  If you get it, try passing in \code{control$abstol} to use delta logLik as the convergence test.}        
      }
  }
  \item{logLik}{Log-likelihood.}
  \item{states}{State estimates from the Kalman filter.}
  \item{states.se}{Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) Shumway & Stoffer.}
  \item{errors}{Any error messages.}
}

\section{Discussion}{
 To ensure that the global maximum-likelihood values are found, it is recommended that initial parameter values be set using Monte Carlo initialization (\code{MLEobj$control$MCInit} = TRUE), particularly if the model is not a good fit to the data.  This requires more compuation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.   For many models, this is unnecessary, but answers should be checked using an initial conditions search before reporting final values.

 \code{MARSSkem()} calls a Kalman filter/smoother (\code{\link{MARSSkf}}) for hidden state estimation.   The algorithm allows two options for the initial state conditions: fixed but unknown or a prior.  In the first case, x0 (whether at t=0 or t=1) is treated as fixed but unknown (estimated); in this case, \code{fixed$V0=0} and x0 is estimated.  This is the default behavior.  In the second case, the initial conditions are specified with a prior and V0!=0.   In the later case, x0 or V0 may be estimated.  MARSS will allow you to try to estimate both, but many researchers have noted that this is not robust so you should fix one or the other.
 
 If you get errors, it generally means that the solution involves an ill-conditioned matrix.    For example, your Q or R matrix is going to a value in which all elements have the same value, for example zero.  If for example, you tried to fit a model with fixed and high R matrix and the variance in that R matrix was much higher than what is actually in the data, then you might drive Q to zero.   Also if you try to fit a structurally inadequate model, then it is not unusual that Q will be driven to zero.  For example, if you fit a model with 1 hidden state trajectory to data that clearly have 2 quite different hidden state trajectories, you might have this problem.  Comparing the likelihood of this model to a model with more structural flexibility should reveal that the structually inflexible model is inadequate (much lower likelihood).
 
     }

\references{ 
  R. H. Shumway and D. S. Stoffer (2006).  Chapter 6 in Time Series Analysis and its Applications.  Springer-Verlag, New York.
  
  Ghahramani, Z. and Hinton, G. E. (1996) Parameter estimation for linear dynamical systems. Technical Report CRG-TR-96-2, University of
Totronto, Dept. of Computer Science. 

Harvey, A. C. (1989) Chapter 5 in Forecasting, structural time series models and the
Kalman filter. Cambridge University Press, Cambridge, UK.
  
The user guide:  Holmes, E. E. and E. J. Ward (2010) Analysis of multivariate time-series
using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   type \code{RShowDoc("UserGuide",package="MARSS")} to see.          

Holmes, E. E. (2010).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models. Link is on the CRAN page under EMDerivation.pdf.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{MARSSmcinit}}, \code{\link{MARSSkf}}, \code{\link{marssMLE}}, \code{\link{MARSSoptim}}
  }  
\examples{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),] 
  #you can use MARSS to construct a proper marssMLE object.
  MLEobj = MARSS(dat, model=list(Q="diagonal and equal", U="equal"), fit=FALSE) 
  #Pass this MLEobj to MARSSkem to do the fit.
  kemfit = MARSSkem(MLEobj)
}


