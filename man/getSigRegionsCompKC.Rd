\name{getSigRegionsCompKC}
\alias{getSigRegionsCompKC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ KCsmart Comparative calculate the signficant regions }
\description{
  Extract the significant regions from a compKC object for a given false discovery rate (FDR).
}
\usage{
getSigRegionsCompKC(compKc, fdr=.01, maxRegionGap=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compKc}{ A compKc object as created by the 'compareSpmCollection' function }
  \item{fdr}{ The false discovery rate to be used to calculate the significantly different regions from the compKc object }
  \item{maxRegionGap}{ The maximum number of sample points that is allowed to fall under the threshold in a continuous significant region}
}
\details{
  The false discovery rate that is set is used to determine the significant regions. When the compKc object was created by the siggenes method the corresponding cutoff is looked up
  in the siggenes results table, otherwise it is calculated from the permuted data.
  The maxRegionGap determines how many sample points can be under this threshold in a continuous significant region.
}

\value{
 Returns a compKcSigRegions object that contains the significant regions for the given FDR in the 'regionTable' slot.
 The method used to determine the cutoff, the fdr and the cutoff itself are stored in their corresponding slots.
 Use 'plot' to visualize the results.
}

\author{ Jorma de Ronde }

\seealso{ 
\code{\link{compareSpmCollection}},  \code{\link{getSigRegionsCompKC}}
}
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spmc1mb <- calcSpmCollection(hsSampleData, hsMirrorLocs, cl=c(rep(0,10),rep(1,10)))
spmcc1mb <- compareSpmCollection(spmc1mb, nperms=3)
spmcc1mbSigRegions <- getSigRegionsCompKC(spmcc1mb)

plot(spmcc1mb, sigRegions=spmcc1mbSigRegions)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
