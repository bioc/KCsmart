\name{compareSpmCollection}
\alias{compareSpmCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ KCsmart Comparative calculate null distribution }
\description{
  Compare the samples of one class in the sample point matrix collection to the samples in the other class and calculate the null distribution
}
\usage{
compareSpmCollection(spmCollection, nperms=20, method=c("siggenes", "perm"), siggenes.args=NULL, altcl=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spmCollection}{ An spmCollection object as created by the 'calcSpmCollection' function}
  \item{nperms}{ The number of permutations to be used to calculate the null distribution }
  \item{altcl}{ Instead of using the class vector from the spmCollection object an alternative vector can be used}
  \item{method}{ The method to be used to calculate the null distribution}
  \item{siggenes.args}{ Optional additional arguments to the siggenes function }
}
\details{
  The method to be used to determine significant regions can either be the SAM methodology from the siggenes package or a signal-to-noise/permutation based method.
  For more information regarding the siggenes method please check the corresponding package. 
}
\value{
 Returns a compKc object which returns the original data and, depending on the method used, the permuted data or the fdr-delta value combinations as calculated by the siggenes package.
}
\references{ }
\author{ Jorma de Ronde }
\note{ 
}
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
