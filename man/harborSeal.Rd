\name{harborSeal}
\alias{harborSeal}
\alias{harborSealnomiss}
\alias{harborSealWA}
\docType{data}
\title{ Harbor Seal Population Count Data (Log counts) }
\description{
 Data sets used in MARSS vignettes in the \code{\link{MARSS-package}}.  These are data sets based on LOGGED count data from Oregon, Washington and California sites where harbor seals were censused while hauled out on land. "harborSeallnomiss" is an extrapolated data set where missing values in the original dataset have been extrapolated so that the data set can be used to demonstrate fitting population models with different underlying structures.
}
\usage{
data(harborSeal)
data(harborSealnomiss)
data(harborSealWA)
}
\format{
  Matrix "harborSeal" contains columns "Years", "StraitJuanDeFuca", "SanJuanIslands", "EasternBays", "PugetSound", "HoodCanal", "CoastalEstuaries", "OlympicPeninsula", "CA.Mainland",  "OR.NorthCoast", "CA.ChannelIslands", and "OR.SouthCoast".
    Matrix "harborSealnomiss" contains columns "Years", "StraitJuanDeFuca", "SanJuanIslands", "EasternBays", "PugetSound", "HoodCanal", "CoastalEstuaries", "OlympicPeninsula", "OR.NorthCoast", and "OR.SouthCoast".
 Matrix "harborSealWA" contains columns "Years", "SJF", "SJI", "EBays", "PSnd", and "HC", representing the same five sites as the first five columns of "harborSeal".
}
\details{
Matrix "harborSealWA" contains the original 1978-1999 LOGGED count data for five inland WA sites. Matrix "harborSealnomiss" contains 1975-2003 data for the same sites as well as four coastal sites, where missing values have been replaced with extrapolated values.  Matrix "harborSeal" contains the original 1975-2003 LOGGED data (with missing values) for the WA and OR sites as well as a CA Mainland and CA ChannelIslands time series.
}
\source{
  Jeffries et al. 2003.  Trends and status of harbor seals in Washington State: 1978-1999. Journal of Wildlife Management 67(1):208-219. 
}
\examples{
str(harborSealWA)
str(harborSealnomiss)
str(harborSeal)
}
\keyword{datasets}
