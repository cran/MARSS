\name{MARSSresiduals}
\alias{MARSSresiduals}

\title{ MARSS Residuals }
\description{
  The normal residuals function is \code{\link[=residuals.marssMLE]{residuals}()}. \code{MARSSresiduals()} returns residuals as a list of matrices while \code{\link[=residuals.marssMLE]{residuals}()} returns the same information in a data frame.  This function calculates the residuals, residuals variance, and standardized residuals for the one-step-ahead (conditioned on data up to \eqn{t-1}), the smoothed (conditioned on all the data), and contemporaneous (conditioned on data up to \eqn{t}) residuals.
}
\usage{
MARSSresiduals(object, ..., type = c("tT", "tt1", "tt"), 
    normalize = FALSE, silent = FALSE, 
    fun.kf = c("MARSSkfas", "MARSSkfss"))
}
\arguments{
  \item{object}{ An object of class \code{marssMLE}.}
  \item{...}{ Additional arguments to be passed to the residuals functions. For type="tT", \code{Harvey=TRUE} can be passed into to use the Harvey et al (1998) algorithm. }
  \item{type}{ \code{"tT"} for smoothed residuals conditioned on all the data \eqn{t=1} to \eqn{T}, aka smoothation residuals. \code{"tt1"} for one-step-ahead residuals, aka innovations residuals. \code{"tt"} for contemporaneous residuals.}
  \item{normalize}{ TRUE/FALSE See details.}
  \item{silent}{ If TRUE, do not print inversion warnings. }
  \item{fun.kf}{ Kalman filter function to use. Can be ignored. }
}
\value{
A list of the following components
  \item{model.residuals}{ The model residuals (data minus model predicted values) as a n x T matrix. }
  \item{state.residuals}{ The  state residuals. This is the state residual for the transition from \eqn{t=t} to \eqn{t+1} thus the last time step will be NA (since \eqn{T+1} is past the data). State residuals do not exist for the \code{type="tt"} case (since this would required the expected value of \eqn{\mathbf{X}_t}{X(t)} conditioned on data to \eqn{t+1}). }
  \item{residuals}{ The residuals as a (n+m) x T matrix with \code{model.residuals} on top and \code{state.residuals} below. }
  \item{var.residuals}{ The variance of the model residuals and state residuals as a (n+m) x (n+m) x T matrix with the model residuals variance in rows/columns 1 to n and state residuals variances in rows/columns n+1 to n+m. The last time step will be all NA since the state residual is for \eqn{t=t } to \eqn{t+1}.}
  \item{std.residuals}{ The Cholesky standardized residuals as a (n+m) x T matrix. This is \code{residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals}. }
  \item{mar.residuals}{ The marginal standardized residuals as a (n+m) x T matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}. }
  \item{bchol.residuals}{ The Block Cholesky standardized residuals as a (n+m) x T matrix. This is \code{model.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[1:n,1:n,]} and \code{state.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[(n+1):(n+m),(n+1):(n+m),]}. }
    \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data. Returned as a n x T matrix.  For observed data, this will be the observed model residuals. For unobserved data, this will be 0 if \eqn{\mathbf{R}}{R} is diagonal but non-zero if \eqn{\mathbf{R}}{R} is non-diagonal. See \code{\link{MARSSresiduals.tT}()}. }
    \item{var.obs.residuals}{ The variance of the model residuals conditioned on the observed data. Returned as a n x n x T matrix.  For observed data, this will be 0. See \code{\link{MARSSresiduals.tT}()}. }
  \item{msg}{ Any warning messages. This will be printed unless Object$control$trace = -1 (suppress all error messages). }

}
\details{

For smoothed residuals, see \code{\link{MARSSresiduals.tT}()}.

For one-step-ahead residuals, see \code{\link{MARSSresiduals.tt1}()}.

For contemporaneous residuals, see \code{\link{MARSSresiduals.tt}()}.

\strong{Standardized residuals}

Standardized residuals have been adjusted by the variance of the residuals at time \eqn{t} such that the variance of the residuals at time \eqn{t} equals 1. Given the normality assumption, this means that one typically sees +/- 2 confidence interval lines on standardized residuals plots. 

\code{std.residuals} are Cholesky standardized residuals. These are the residuals multiplied by the inverse of the lower triangle of the Cholesky decomposition of the variance matrix of the residuals: 
\deqn{ \hat{\Sigma}_t^{-1/2} \hat{\mathbf{v}}_t.}{ hatSigma_t^{-1/2} hatv_t.}
These residuals are uncorrelated with each other, although they are not necessarily temporally uncorrelated (innovations residuals are temporally uncorrelated). 

The interpretation of the Cholesky standardized residuals is not straight-forward when the \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in \eqn{\textrm{MVN}(0,\mathbf{I})}{MVN(0,I)} space.  For example, if v is 2x2 correlated errors with variance-covariance matrix \eqn{\mathbf{R}}{R}. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the Cholesky standardized residuals can look rather non-intuitive.

\code{mar.residuals} are the marginal standardized residuals. These are the residuals multiplied by the inverse of the diagonal matrix formed from the square-root of the diagonal of the variance matrix of the residuals: 
\deqn{ \textrm{dg}(\hat{\Sigma}_t)^{-1/2} \hat{\mathbf{v}}_t,}{ dg(hatSigma_t)^{-1/2} hatv_t,} where \eqn{dg(A)} is the square matrix formed from the diagonal of \eqn{A}, aka \code{diag(diag(A))}. These residuals will be correlated if the variance matrix is non-diagonal. 

The Block Cholesky standardized residuals are like the Cholesky standardized residuals except that the full variance-covariance matrix is not used, only the variance-covariance matrix for the model or state residuals (respectively) is used for standardization. For the model residuals, the Block Cholesky standardized residuals will be the same as the Cholesky standardized residuals because the upper triangle of the lower triangle of the Cholesky decomposition (which is what we standardize by) is all zero. For \code{type="tt1"} and \code{type="tt"}, the Block Cholesky standardized state residuals will be the same as the Cholesky standardized state residuals because in the former, the model and state residuals are uncorrelated and in the latter, the state residuals do not exist. For \code{type="tT"}, the model and state residuals are correlated and the Block Cholesky standardized residuals will be different than the Cholesky standardized residuals.

\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} are \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q} and the model is assumed to be written as
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{v}_t}{ y(t) = Z x(t) + a + v(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{w}_t}{ x(t) = B x(t-1) + u + w(t)}
If \code{normalize=TRUE}, the model is assumed to be written as
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{H}\mathbf{v}_t}{ y(t) = Z x(t) + a + Hv(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{G}\mathbf{w}_t}{ x(t) = B x(t-1) + u + Gw(t)}
with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} equal to \eqn{\mathbf{I}}{I} (identity).

\strong{Missing or left-out data}

See the discussion of residuals for missing and left-out data in \code{\link{MARSSresiduals.tT}()}.

}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
} 
\seealso{ \code{\link{residuals.marssMLE}()}, \code{\link{MARSSresiduals.tT}()}, \code{\link{MARSSresiduals.tt1}()}, \code{\link{plot.marssMLE}()} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  fit <- MARSS(dat)
  
  #state smoothed residuals
  state.resids1 <- MARSSresiduals(fit, type="tT")$state.residuals
  #this is the same as
  states <- fit$states
  Q <- coef(fit, type="matrix")$Q
  state.resids2 <- states[,2:30]-states[,1:29]-matrix(coef(fit,type="matrix")$U,2,29)
  #compare the two
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #normalize to variance of 1
  state.resids1 <- MARSSresiduals(fit, type="tT", normalize=TRUE)$state.residuals
  state.resids2 <- (solve(t(chol(Q))) \%*\% state.resids2)
  cbind(t(state.resids1[,-30]), t(state.resids2))

  #one-step-ahead standardized residuals
  MARSSresiduals(fit, type="tt1")$std.residuals
}
\references{
Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045.

See also the discussion and references in \code{\link{MARSSresiduals.tT}()},  \code{\link{MARSSresiduals.tt1}()} and \code{\link{MARSSresiduals.tt}()}.
}

