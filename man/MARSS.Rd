\name{MARSS}
\alias{MARSS}
\alias{MARSSsettings}
\title{ Interface MARSS Model Specification and Estimation }
\description{
  A top-level \code{\link{MARSS-package}} function to perform model specification and estimation for multivariate autoregressive state-space (MARSS) models.  To open the user guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}.  To open an overview page with package information, type \code{RShowDoc("index",package="MARSS")} or go to \code{\link{MARSS-package}}.    

The MARSS package fits state-space models that can be written in the \code{\link{marssm}} form:
\describe{
\item{x(t+1) = (t(x(t)) (x) I_m)(f_b(t)+D_b(t)beta) + (f_u(t)+D_u(t)upsilon) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = (t(x(t)) (x) I_n)(f_z(t)+D_z(t)zeta) + (f_z(t)+D_z(t)alpha) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}{} 
}
where beta, upsilon, zeta, and alpha are column vectors of estimated values, the f are column vectors of inputs, and the D are matrices of inputs.  The f and D are potentially time-varying. (x) means kronecker product and I_p is a p x p identity matrix.  The function \code{MARSS()} is used to fit MARSS models using the argument \code{form} to specify the type of state-space model being fit.  \code{MARSS()} reformulates the model into marssm form and fits the model.

The default MARSS form is "marxss" which is the state-space model:
\describe{
\item{x(t) = B(t) x(t-1) + u(t) + C(t)c(t) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + a(t) + D(t)d(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}{} 
}
If you are working with models with time-varying parameters, it is important to notice the time-index for the parameters in the process equation (the x equation).  In some formulations (e.g. in the \code{\link[KFAS]{KFAS}}), the process equation is x(t)=B(t-1)x(t-1)+w(t-1) so B(t-1) goes with x(t) not B(t). Thus one needs to be careful to line up the time indices when passing in time-varying parameters to MARSS.
}
\usage{
MARSS(y, 
    inits=NULL,
    model=NULL,
    miss.value=NA,
    method = "kem",
    form = "marxss",
    fit=TRUE, 
    silent = FALSE,
    control = NULL,
    MCbounds = NULL,
    fun.kf = "MARSSkfas",
    ...) 
}

\arguments{
The default settings for the optional arguments are set in \code{MARSSsettings.R} and are given below in the details section.  
  \item{y}{ A n x T matrix of n time series over T time steps. }
  \item{inits}{ List with up to 8 matrices Z, A, R, B, U, Q, x0, V0, specifying initial values for beta, upsilon, zeta and alpha. All matrices have 1 column.  Each row is one estimtated value.  This is the \code{$par} element of a \code{\link{marssMLE}} object and is in \code{\link{marssm}} form (\code{form="marss"}).  Note the output of \code{coef(marssMLE,type="list")} is similar but is the parameter list in the form specified by \code{marssMLE$form}, which is probably not \code{marss}.
    \itemize{
    \item \code{B }  beta: Initial value(s) for estimates values in the B matrix.  vec(B(1))=f_b(1)+D_b(1)beta_init
    \item \code{U }  upsilon: Initial value(s) for U matrix. U(1)=f_u(1)+D_u(1)upsilon_init 
    \item \code{Q }  Initial value(s) for process error variance(s). vec(Q(1))=f_q(1)+D_q(1)q_init 
    \item \code{Z }  zeta: Initial value(s) for estimates values in the Z matrix.  vec(Z(1))=f_z(1)+D_z(1)zeta_init
    \item \code{A }  alpha: Initial value(s) for observation bias. A(1)=f_a(1)+D_a(1)alpha_init 
    \item \code{R }  Initial value(s) for non-process (observation) error variance(s). vec(R(1))=f_r(1)+D_r(1)r_init
    \item \code{x0} Initial value(s) for hidden state(s) at time=1 or 0.  
    \item \code{V0} Initial variance(s) for hidden state(s) at time=1 or 0. 
    }
  }
  \item{model}{ Model specification using parameter model descriptions. See Details. }
  \item{miss.value}{ How are missing values represented in the data? }
  \item{method}{ Estimation method. MARSS allows method="kem" (see \code{\link{MARSSkem}} for more info) and "BFGS" (see \code{\link{MARSSoptim}} for more info). }
  \item{form}{ The equation form used to translate the model argument into a marssm model object for the fitting algorithms.  The default is \code{form="marxss"} (see \code{\link{MARSS.marxss}}).  \code{form="dfa"} is also available for constructing Dynamic Factor Analysis (DFA) models (see \code{\link{MARSS.dfa}}). }
  \item{fit}{ TRUE/FALSE Whether to fit the model to the data. If FALSE, a marssMLE object with only the model is returned. }
  \item{silent}{ TRUE/FALSE Suppresses printing of full error messages, warnings, progress bars and convergence information. }
  \item{control}{ Estimation options for the maximization algorithm. The typically used control options for method="kem" are below but see  \code{\link{marssMLE}} for the full list of control options.  Note many of these are not allowed if method="BFGS"; see  \code{\link{MARSSoptim}} for the allowed control options for this method.
    \itemize{
      \item{\code{minit}}{   The minimum number of iterations to do in the maximization routine (if needed by method). If \code{method="kem"}, this is an easy way to up the iterations and see how your estimates are converging. (positive integer)} 
      \item{\code{maxit}}{   Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). }
      \item{\code{min.iter.conv.test}}{ Minimum iterations to run before testing convergence via the slope of the log parameter versus log iterations.} 
      \item{\code{conv.test.deltaT=9}}{ Number of iterations to use for the testing convergence via the slope of the log parameter versus log iterations. }
      \item{\code{conv.test.slope.tol}}{ The slope of the log parameter versus log iteration to use as the cut-off for convergence.  The default is 0.5 which is a bit high.  For final analyses, this should be set lower. }
      \item{\code{abstol}}{ Optional convergence tolerance for the maximization routine. If set to something other than NULL, then the absolute tolerance will be used as the convergence test.  Otherwise, the slope of the log of the parameters versus the log iteration is used as the test.} 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{An integer specifying the level of information recorded and error-checking run during the algorithms.  trace=0, specifies basic error-checking and brief error-messages; trace>0 will print full error messages. In addition if trace>0, the Kalman filter output will be added to the outputted marssMLE object.  Additional information recorded depends on the method of maximization. For the EM algorithm, a record of each parameter estimate for each EM iteration will be added.   See \code{\link{optim}} for trace output details for the BFGS method.   trace=-1 will turn off most internal error-checking and most error messages.  The internal error checks are time expensive so this can speed up MARSS.  This is particularly useful for bootstrapping and simulation studies.}
      \item{\code{MCInit}}{   If TRUE, do a Monte Carlo search of the initial condition space. (T/F) } 
      \item{\code{numInits}}{   Number of random initial value draws if MCInit=TRUE (ignored otherwise). (positive integer)} 
      \item{\code{numInitSteps}}{   Number of EM iterations for each random initial value draw if MCInit=TRUE (ignored otherwise). (positive integer)}
      \item{\code{silent}}{   1 or TRUE, Suppresses all printing including progress bars, error messages and convergence information. 0, Turns on all printing of progress bars, fitting information and error messages.  2, Prints a brief success/failure message.    }
  } }
    \item{MCbounds}{   Bounds on the uniform distributions from which initial values will be drawn if MCInit=TRUE (ignored otherwise). The Q (and analogously R) random matrices are created by a random draw from a wishart distribution where df=bounds[1] and S=diag(bounds[2],m). }
    \item{fun.kf}{   What Kalman filter function to use.  MARSS has two: \code{\link{MARSSkfas}} which is based on the Kalman filter in the KFAS package based on Koopman and Durbin and \code{\link{MARSSkfss}} which is a native R implementation of the Kalman filter and smoother in Shumway and Stoffer. The KFAS filter is much faster.  \code{\link{MARSSkfas}} modifies the input and output in order to output the lag-one covariance smoother needed for the EM algorithm (per page 321 in Shumway and Stoffer (2000).}
    \item{...}{ Optional arguments passed to function specified by form. }
}
\details{
 MARSS provides an interface to the base \code{\link{MARSS-package}} functions and allows specification and fitting of MARSS models. The available estimation methods are maximum-likelihood via an EM algorithm (\code{method="kem"}) or via a quasi-Newton algorithm provided by function \code{optim} (\code{method="BFGS"}).  The function \code{MARSS()} allows the user to specify models using the \code{model} argument which is a list specifying the structure for each parameter matrix.  The user can use text strings for common classes of parameter matrices.  It also allows the user to specify fixed values for matrices by passing in numeric matrices in the \code{model} list.  Models with a combination of fixed values and shared estimated values can be specified using list matrices. See \code{\link{marssm}} or the user guide (reference and link below) for documentation and instructions on specifying list matrices for the model parameters.
 
 Many different types of multivariate time-series models can be converted to the MARSS form (see the User Guide or Quick Start guide). \code{MARSS()} also the user to specify the form of their model using the argument \code{form}. The default form is a multivariate lag-1 (time-varying) state-space model which is denoted using \code{form="marxss"}. Look at \code{\link{MARSS.marxss}}) to see how the model argument for the marxss form is specified.   The model argument is a list, but the elements of that and how it is converted to a marssm model object depends intimately on the equation form.  Thus you will need to refer to that appropriate help page (\code{MARSS.formname}) for your equation form.  

 The likelihood surface for MARSS models can be  multimodal or with strong ridges.  It is  recommended that for final analyses the ML estimates are checked by using the Monte Carlo initial conditions search using \code{MCInit=TRUE} in the \code{control} list.  This requires more computation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.  Also it is wise to check the EM results against the BFGS results (if possible) since if there are strong ridges in the likelihood, this seems to slow down the EM algorithm considerably and cause the algorithm to report convergence far from the ML values.  EM steps up the likelihood and the convergence test is based on the rate of change of the LL in each step; once on a strong ridge, the steps can slow dramatically.  You can force the algorithm to keep working by setting \code{minit}.  BFGS seems less hindered by the ridges but can be prodigiously slow for some multivariate problems.
 
Note: One of the most common inquiries we get is "Why doesn't my model converge?" or "Why is Q (or R) going to 0?".  There are two main reasons for this. 1) The model is structurally indeterminate. This means that there are multiple equivalent solutions.  If you were to look at the likelihood surface, you would see ridges. 2) The model is structurally inadequate. This means you specified a model that cannot (naturally) fit the data and since you told MARSS to fit that data with a bad model, it does the best it can.  Some typical inadequacies we see are a) the model implies that the data have a mean of 0, yet the data do not have a mean of zero, 2) The data are stationary, yet the model would not natually generate a stationary time series (because the state model is non-stationary), and 3) the data are not-stationary yet the model implies stationary data (and would not natually generate an unstationary time series).
}

\value{
  An object of class marssMLE with the following components:
  \item{model}{ MARSS model specification. }
  \item{start}{ List with 8 matrices A, R, B, U, Q, Z, x0, V0, specifying initial values that were used for the values estimated. } 
  \item{control}{ A list of estimation options, as specified by arguments \code{control}. }
  \item{method}{ Estimation method. }
  If \code{fit=TRUE}, the following are also added to the marssMLE object.
  If \code{fit=FALSE}, an marssMLE object ready for fitting via the specified \code{method} is returned.
  \item{par}{ A list with 8 matrices of estimated parameter values Z, A, R, B, U, Q, x0, V0. The model argument is needed to translate this to parameter matrices as it is a vector of only the estimated elements in a model parameter. }
  \item{kf}{ A list containing Kalman filter/smoother output from \code{\link{MARSSkf}}.}
  \item{numIter}{ Number of iterations required for convergence. }
  \item{convergence}{ Convergence status. 0 means converged successfully.  Anything else is a warning or error. 2 means the MLEobj has an error; the MLEobj is returned so you can debug it.  The other numbers are errors during fitting.  The error code depends on the fitting method.  See \code{\link{MARSSkem}} and \code{\link{MARSSoptim}}. }
  \item{logLik}{ Log-likelihood. }
  \item{AIC}{ Akaike's Information Criterion. }
  \item{AICc}{ Sample size corrected AIC. }
}
\references{ 
The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.          

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models.  Type \code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.
}
\author{ 
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, kellie(dot)wills(at)noaa(dot)gov
}
\seealso{
\code{\link{marssm}} \code{\link{marssMLE}} \code{\link{MARSSkem}} \code{\link{MARSSoptim}} \code{\link{MARSS-package}}
}
\examples{
#harborSealWA is a n=5 matrix of logged population counts
dat = t(harborSealWA)
dat = dat[2:4,] #remove the year row
#fit a model with 1 hidden state and 3 observation time series
kemfit = MARSS(dat, model=list(Z=factor(c(1,1,1)), 
   R="diagonal and equal"))
kemfit$model #This gives a description of the model
print(kemfit$model) # same as kemfit$model
summary(kemfit$model) #This shows the model structure

#add CIs to a marssMLE object
#default uses an estimated Hessian matrix
kem.with.hess.CIs = MARSSparamCIs(kemfit)   
kem.with.hess.CIs #print with se's and CIs

#fit a model with 3 hidden states (default)
kemfit = MARSS(dat, silent=TRUE)  #suppress printing
kemfit #print information on the marssMLE object

#fit a model with 3 correlated hidden states 
# with one variance and one  covariance
#maxit set low to speed up example, but more iterations are need for convergence  
kemfit = MARSS(dat, model=list(Q="equalvarcov"), control=list(maxit=50))

#fit a model with 3 correlated hidden states
#maxit set low to speed up example, but more iterations are need for convergence  
kemfit = MARSS(dat, model=list(Q="unconstrained"), control=list(maxit=50))

#fit a model with 3 independent hidden states 
#where each observation time series is independent
#the hidden trajectories 1-3 & 4-5 share their U parameter
kemfit = MARSS(dat, model=list(U=matrix(c("N","S","S"),3,1)))

#same model, but with fixed independent observation errors
#and the 3rd x processes are forced to have a U=0
#Notice how a list matrix is used to combine fixed and estimated elements
#all parameters can be specified in this way using list matrices
kemfit = MARSS(dat, model=list(U=matrix(list("N","N",0),3,1), R=diag(0.01,3)))

#fit a model with 2 hidden states (north and south) 
#where observation time series 1-2 are north and 3 is south
#Make the hidden state process independent with same process var
#Make the observation errors different but independent
#Make the growth parameters (U) the same
kemfit = MARSS(dat, model=list(Z=factor(c("N","N","S")),
   Q="diagonal and equal",R="diagonal and unequal",U="equal"))

#print the model followed by the marssMLE object
kemfit$model
kemfit

#simulate some new data from our fitted model
sim.data=MARSSsimulate(kemfit, nsim=10, tSteps=100)


\dontrun{
#Compute bootstrap AIC for the model; this takes a long, long time
kemfit.with.AICb = MARSSaic(kemfit, output = "AICbp")
kemfit.with.AICb
}
}

