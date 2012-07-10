\name{plankton}
\alias{ivesDataLP}
\alias{ivesDataByWeek}
\alias{lakeWAplankton}
\docType{data}
\title{ Plankton Data Sets }
\description{
  Example data sets for use in MARSS vignettes for the \code{\link{MARSS-package}}. 
  
  The "lakeWAplankton" plankton counts have been standardized to a mean of zero and variance of 1 (Z-score transformation). The second column in "lakeWAplankton" is simply a index for the count.  Columns 3 and 4 are month and month^2 but they have also been Z-score transformed.  Since MARSS functions require time to be across columns, these data matrices must be transposed before passing into MARSS functions.
  
  The Ives data are from Ives et al. (2003) for West Long Lake (the low planktivory case).  The Ives data are unlogged.  \code{ivesDataLP} and \code{ivesDataByWeek} are the same data with LP having the missing weeks in winter removed while in ByWeek, the missing values are left in.  The phosporous column is the experimental input rate + the natural input rate for phosphorous, and Ives et al. used 0.1 for the natural input rate when no extra phosporous was added.  The phosporous input rates for weeks with no sampling (and no experimental phosporous input) have been filled with 0.1 in the "by week" data.
}
\usage{
data(ivesDataLP)
data(ivesDataByWeek)
data(lakeWAplankton)
}
\format{
The data are provided as a matrix with time running down the rows.
}
\source{
\itemize{
\item{ivesDataLP and ivesDataByWeek}{Ives, A. R. Dennis, B. Cottingham, K. L. Carpenter, S. R. (2003) Estimating community stability and ecological interactions from time-series data. Ecological Monographs, 73, 301-330. }
\item{lakeWAplankton}{Hampton, S. E. Scheuerell, M. D. Schindler, D. E. (2006)  Coalescence in the Lake Washington story: Interaction strengths in a planktonic food web. Limnology and Oceanography, 51, 2042-2051. }
}
}
\examples{
str(ivesDataLP)
str(ivesDataByWeek)
str(lakeWAplankton)
}
\keyword{datasets}
