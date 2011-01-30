\name{graywhales}
\alias{graywhales}
\alias{grouse}
\alias{isleRoyal}
\alias{prairiechicken}
\alias{wilddogs}
\alias{rockfish}
\alias{okanaganRedds}
\alias{kestrel}
\docType{data}
\title{ Population Data Sets }
\description{
  Example data sets for use in MARSS PVA vignettes in the \code{\link{MARSS-package}} user guide. The data sets are matrices with year in the first column and counts in other columns.  Since MARSS functions require time to be across columns, these data matrices must be transposed before passing into MARSS functions.
}
\usage{
data(graywhales)
data(grouse)
data(isleRoyal)
data(prairiechicken)
data(wilddogs)
data(kestrel)
data(okanaganRedds)
data(rockfish)
}
\format{
The data are supplied as a matrix with years in the first column and counts in the second (and third for isleRoyal) columns.
}
\source{
\itemize{
\item{graywhales }{ Gerber L. R., Master D. P. D. and Kareiva P. M. (1999) Gray whales and the value of monitoring data in implementing the U.S. Endangered Species Act.  Conservation Biology, 13, 1215-1219. UNLOGGED counts. }
\item{grouse }{ Hays D. W., Tirhi M. J. and Stinson D. W. (1998) Washington state status report for the sharptailed grouse. Washington Department Fish and Wildlife, Olympia, WA. 57 pp. UNLOGGED counts. }
\item{isleRoyal }{ Peterson R. O., Thomas N. J., Thurber J. M., Vucetich J. A. and Waite T. A. (1998) Population limitation and the wolves of Isle Royale. In: Biology and Conservation of Wild Canids (eds. D. Macdonald and C. Sillero-Zubiri). Oxford University Press, Oxford, pp. 281-292. UNLOGGED counts. }
\item{prairiechicken }{ Peterson M. J. and Silvy N. J. (1996) Reproductive stages limiting productivity of the endangered Attwater's prairie chicken. Conservation Biology, 10, 1264-1276. UNLOGGED counts. }
\item{wilddogs }{ Ginsberg, J. R., Mace, G. M. and Albon, S. (1995). Local extinction in a small and declining population: Wild Dogs in the Serengeti. Proc. R. Soc. Lond. B, 262, 221-228. UNLOGGED population counts. }
\item{okanaganRedds }{ A dataset of Chinook salmon redd (egg nest) surveys.  This data comes from the Okanagan River in Washington state, a major tributary of the Columbia River (headwaters in British Columbia). UNLOGGED redd counts. }
\item{rockfish }{ LOGGED catch per unit effort data for Puget Sound total total rockfish (mix of species) from a series of different types of surveys.}
\item{kestrel }{ Three time series of American kestrel logged abundance from adjacent Canadian provinces along a longitudinal gradient (British Columbia, Alberta, Saskatchewan).  Data have been collected annually, and corrected for changes in observer coverage and detectability. LOGGED abundance estimates.}
}
}
\examples{
str(graywhales)
str(grouse)
str(isleRoyal)
str(prairiechicken)
str(wilddogs)
str(kestrel)
str(okanaganRedds)
str(rockfish)
}
\keyword{datasets}
