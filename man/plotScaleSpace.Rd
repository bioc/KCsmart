\name{plotScaleSpace}
\alias{plotScaleSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot multiple significant regions in one figure }
\description{
  Plots significant regions in different scale spaces in one figure
}
\usage{
plotScaleSpace(spms, sigLevels, chromosomes=NULL, type='b')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spms}{ List of sample point matrices }
  \item{sigLevels}{ List of significance levels }
  \item{chromosomes}{ Takes a vector of chromosomes to be plotted. Defaults to all chromosomes. }
  \item{type}{ Determines which data is plotted. 'g' for gains only, 'l' for losses only and 'b' for both. When type='b' is used, two devices (x11) will be opened.}
}
\details{
  Takes sample point matrices that were calculated using (different) kernel widths (sigma), then calculates the significant regions given the cutoffs as defined by 'sigLevels' and plots these in one figure.
}
\value{
  Depending on the 'type' parameter, produces one or two plots, one for the gains and one for the losses. The heatmap color indicates the level of the gain or loss.
}
\author{ Jorma de Ronde }
\note{ 
}
\seealso{ \code{\link{plot}} }
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)
spm4mb <- calcSpm(hsSampleData, hsMirrorLocs, sigma=4000000)

siglevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)
siglevel4mb <- findSigLevelTrad(hsSampleData, spm4mb, n=3)

plotScaleSpace(list(spm1mb, spm4mb), list(siglevel1mb, siglevel4mb), type='g')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
