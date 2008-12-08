\name{plot}
\alias{plot}
\alias{plot,scaleSpace,missing-method}
\alias{plot,samplePointMatrix,missing-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a sample point matrix }
\description{
Plot the sample point matrix or parts of it
}
\usage{
plot(x,y, ...)
\S4method{plot}{scaleSpace,missing}(x, y, spm, type='b', ...)
\S4method{plot}{samplePointMatrix,missing}(x, y, type="b", sigLevels=NULL, chromosomes=NULL, colinAxis=NULL, fillColor=NULL, maploc=NULL, interpolation=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{either an object of class samplePointMatrix or scaleSpace}
  \item{y}{object of class missing}
  \item{type}{ Determines which data is plotted. 'g' for gains only, 'l' for losses only and 'b' and '1' for both in one plot device}
  \item{spm}{ add stuff here}
  \item{sigLevels}{ If given, the cutoffs will be drawn as lines in the plots. Optional }
  \item{chromosomes}{ Takes a vector of chromosomes to be plotted. Defaults to all chromosomes. }
  \item{colinAxis}{ Allows you to override default behaviour of axis labeling. Choose False for genomic position labeling for each individual chromosome, True for colinear labeling.}
  \item{fillColor}{ Allows you to choose the colors used to fill the significant areas under the curve. Takes a list with the 'pos' element giving the color for the gains and the 'neg' element the color for the losses. }
  \item{maploc}{ Currently not in use}
  \item{interpolation}{ Determines which points from the sample point matrix will actually be plotted. If the value of 'interpolation' is n, then every n-th point will be plotted. The default value of 1 will results in all points being plotted. This can be useful when a high density sample point matrix results in big file size when exporting the image (especially to pdf or eps format).}
  \item{...}{Any other parameters you would like to pass to 'plot'. See 'par' for more details.}
}
\details{
}
\value{
 Plots the sample point matrix. The gains and the losses are plotted separately. The KC normalized score is plotted on the y-axis, the genomic position on the x-axis.
 If centromeres are present these are represented by dotted, lightblue lines.
 Setting type to 'b' or to '1' will both make the plot appear in one plot device, '1' will plot the gains and the losses in one plot, 'b' will plot the gains and losses separately.
The function 'idPoints' can be used to identify points in the sample point matrix plot. See the corresponding documentation for details.
}
\author{ Jorma de Ronde }
\note{ 
}
\seealso{  \code{\link{calcSpm}}, \code{\link{plotScaleSpace}}, \code{\link{idPoints}} }
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)

plot(spm1mb)
plot(spm1mb, interpolation=10)
plot(spm1mb, chromosomes=c(1,4,'X'))

siglevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)
plot(spm1mb, chromosomes=c(1,4,'X'), sigLevels=siglevel1mb)
plot(spm1mb, chromosomes=c(1,4,'X'), sigLevels=siglevel1mb, fillColor=list(pos='darkred',neg='darkgreen'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}