\name{MARSSaic}
\alias{MARSSaic}
\title{ AIC for MARSS models }
\description{
  Calculates AIC, AICc, a parametric bootstrap AIC (AICbp) and a non-parametric bootstrap AIC (AICbb).  This is a base function in the \code{\link{MARSS-package}}.
}
\usage{
MARSSaic(MLEobj, output = c("AIC", "AICc"), 
  Options = list(nboot = 1000, return.logL.star = FALSE, 
  silent = FALSE))
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{MARSSkem()}. }
  \item{output}{ A vector containing one or more of the following: "AIC", "AICc", "AICbp", "AICbb", "AICi", "boot.params". See Details.}
  \item{Options}{ A list containing:
    \itemize{
      \item \code{nboot}  Number of bootstraps (positive integer) 
      \item \code{return.logL.star}  Return the log-likelihoods for each bootstrap? (T/F)
      \item \code{silent}  Suppress printing of the progress bar during AIC bootstraps? (T/F)  
    } 
  }
}
\details{
  When sample size is small, Akaike's Information Criterion (AIC) under-penalizes more complex models. The most commonly used small sample size corrector is AICc, which uses a penalty term of Kn/(n-K-1), where K is the number of estimated parameters. However, for time series models, AICc still under-penalizes complex models; this is especially true for MARSS models.

Two small-sample estimators specific for MARSS models have been developed. Cavanaugh and Shumway (1997) developed a variant of bootstrapped AIC using Stoffer and Wall's (1991) bootstrap algorithm ("AICbb"). Holmes and Ward (2010) developed a variant on AICb ("AICbp") using a parametric bootstrap. The parametric bootstrap permits AICb calculation when there are missing values in the data, which Cavanaugh and Shumway's algorithm does not allow. More recently, Bengtsson and Cavanaugh (2006) developed another small-sample AIC estimator, AICi, based on fitting candidate models to multivariate white noise.

When \code{output} includes both \code{"AICbp"} and \code{"boot.params"}, the bootstrapped parameters from \code{"AICbp"} will be added to \code{MLEobj}.  
}
\value{
  Returns the \code{\link{marssMLE}} object that was passed in with additional AIC components added on top as specified in the 'output' argument.
}
\references{
Holmes, E. E. and E. J. Ward. 2010. Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112; this is the user manual accesses via \code{show.doc(MARSS, manual)}
 
Bengtsson, T., and J. E. Cavanaugh. 2006. An improved Akaike information criterion for state-space model selection. Computational Statistics & Data Analysis 50:2635-2654.

Cavanaugh, J. E., and R. H. Shumway. 1997. A bootstrap variant of AIC for state-space model selection. Statistica Sinica 7:473-496.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ \code{\link{MARSSboot}} }
\examples{
  dat = t(harborSealWA)
  dat = dat[2:nrow(dat),]
  kem = MARSS(dat, constraint=list(Z=factor(c(1,1,1,1,1)), 
    R="diagonal and unequal"), control=list(minit=100))
  kemAIC = MARSSaic(kem, output=c("AIC","AICc"))
}

