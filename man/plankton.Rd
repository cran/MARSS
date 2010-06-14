\name{plankton}
\alias{ivesDataLP}
\alias{lakeWAplankton}
\docType{data}
\title{ Plankton Data Sets }
\description{
  Example data sets for use in MARSS vignettes for the \code{\link{MARSS-package}}. The "lakeWAplankton" plankton counts have been standardized to a mean of zero and variance of 1 (Z-score transformation). The second column in "lakeWAplankton" is simply a index for the count.  Columns 3 and 4 are month and month^2 but they have also been Z-score transformed.  Since MARSS functions require time to be across columns, these data matrices must be transposed before passing into MARSS functions.
}
\usage{
data(ivesDataLP)
data(lakeWAplankton)
}
\format{
The data are provided as a matrix with time running down the rows.
}
\source{
\itemize{
\item{ivesDataLP}{Ives, A. R. Dennis, B. Cottingham, K. L. Carpenter, S. R. (2003) Estimating community stability and ecological interactions from time-series data. Ecological Monographs, 73, 301-330. }
\item{lakeWAplankton}{Hampton, S. E. Scheuerell, M. D. Schindler, D. E. (2006)  Coalescence in the Lake Washington story: Interaction strengths in a planktonic food web. Limnology and Oceanography, 51, 2042-2051. }
}
}
\examples{
str(ivesDataLP)
str(lakeWAplankton)
}
\keyword{datasets}
