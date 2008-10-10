\name{findSigLevelFdr}
\alias{findSigLevelFdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ This function has not been properly implemented yet }
\description{
  Method to find the cutoff at which gains and losses are considered significant using permutations
}
\usage{
findSigLevelFdr(data, observedSpm, n = 1, fdrTarget=0.05, maxmem=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ aCGH data in the same format as used for 'calcSpm' }
  \item{observedSpm}{ A sample point matrix as produced by 'calcSpm'}
  \item{n}{ Number of permutations }
  \item{fdrTarget}{ Target False Discovery Rate (FDR) }
  \item{maxmem}{ This parameter controls memory usage, set to lower value to lower memory consumption }
}
\details{
  The number of permutations needed for reliable results depends on the data and can not be determined beforehand. As a general rule-of-thumb around 100 permutations should be used for 'quick checks'
  and around 2000 permutations for more rigorous testing. 
  The FDR method is less conservatie than the p-value based approach since instead of controlling the family wise error rate (FWER, P(false positive > 1)) it controls the false
  discovery rate (FDR) (false positives / total number of called data points).
}
\value{
  A list with the cutoffs corresponding to the given FDR

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
