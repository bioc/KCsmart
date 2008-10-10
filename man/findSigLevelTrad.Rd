\name{findSigLevelTrad}
\alias{findSigLevelTrad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find significance level }
\description{
  Method to find the cutoff at which gains and losses are considered significant using permutations
}
\usage{
findSigLevelTrad(data, observedSpm, n = 1, p = 0.05, maxmem = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ aCGH data in the same format as used for 'calcSpm' }
  \item{observedSpm}{ A sample point matrix as produced by 'calcSpm'}
  \item{n}{ Number of permutations }
  \item{p}{ Alpha level for significance }
  \item{maxmem}{ This parameter controls memory usage, set to lower value to lower memory consumption }
}
\details{
  The number of permutations needed for reliable results depends on the data and can not be determined beforehand. As a general rule-of-thumb around 100 permutations should be used for 'quick checks'
  and around 2000 permutations for more rigorous testing.
  
  p is the uncorrected alpha level, the method corrects for multiple testing internally using simple Bonferroni correction. See the referenced publication for more details.
}
\value{
  A list with the cutoffs corresponding to the given alpha level

  \item{pos }{The cutoff for the gains}
  \item{neg }{The cutoff for the losses'}
}
\author{ Jorma de Ronde }

\seealso{\code{\link{plotScaleSpace}} }
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)

sigLevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)

plot(spm1mb, sigLevels=sigLevel1mb)
plotScaleSpace(list(spm1mb), list(sigLevel1mb), type='g')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
