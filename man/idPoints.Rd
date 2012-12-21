\name{idPoints}
\alias{idPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identify points in sample point matrix plot }
\description{
   Identify points in sample point matrix plot 
}
\usage{
idPoints(spm, mode='pos', dev=2, chromosomes=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ The sample point matrix object of which points are to be identified }
  \item{mode}{ Determines which points will be identified: mode='pos' will identify points in gained regions, mode='neg' will identify points in lost regions }
  \item{dev}{ The device on which the sample point matrix was plotted}
  \item{chromosomes}{ If not all chromosomes contained in the sample point matrix were plotted (using the 'chromosomes' argument in the 'plot' command), the same chromosomes must be entered here as an argument}
}
\details{
   Using the mouse pointer points in a sample point matrix plot can be identified by left-clicking on the to-be-identified points. Right-clicking exits the selection and returns the selected points.
}
\value{
  Returns a data.frame listing the the position and the KC score for each identified point.
    \item{KCscore}{ KCscore of the identified point}
    \item{chromosome}{Chromosome on which the identified point is located}
    \item{chromPosition}{Position on the chromosome of the identified point}
    \item{colin}{Co-linear location of the identified point (given the selected chromosomes)}
}
\author{ Jorma de Ronde }

\seealso{ \code{\link{plot}} }
\examples{
data(hsSampleData)
data(hsMirrorLocs)

#spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)

#plot(spm1mb, type=1)
#idPoints(spm1mb)

#x11()
#plot(spm1mb, chromosomes=c(1,2,5))
#idPoints(spm1mb, mode='neg', dev=3, chromosomes=c(1,2,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{iplot}
