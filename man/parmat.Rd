\name{parmat}
\alias{parmat}
\title{ Retrieve Parameter Matrix }
\description{
  Function to return the parmater matrix with fixed and free elements given a model object (\code{\link{is.marssm}}) and the \code{par} element of a marss MLE object (\code{\link{is.marssMLE}}).
}
\usage{
parmat(MLEobj, elem=c("B","U","Q","Z","A","R","x0","V0"), t=1, dims=NULL)
}
\arguments{

  \item{MLEobj}{ A marss MLE object }
  \item{elem}{ The name (character string) of a parameter, e.g. "A". }
  \item{t}{ The time index or third dimension of a 3D matrix }
  \item{dims}{ Tells the function the dimension of elem.  Must be a list with elem as the element name. }
}
\details{
  \code{parmat} constructs a 2D parameter matrix from a marssm model object and \code{par} element of a marss MLE object.  The equation for the vec-ed version of the matrix is \code{marssm$fixed[[elem]]} + \code{marssm$free[[elem]]} x \code{marssMLE$par[[elem]]}.  This is then unvec-ed and the matrix returned for element \code{elem}. 
}
\value{
  The parameter matrix for element \code{elem}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}

