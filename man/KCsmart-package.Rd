\name{KCsmart-package}
\alias{KCsmart-package}
\alias{KCsmart}
\docType{package}
\title{
KCsmart
}
\description{
Multiple sample aCGH analysis using kernel convolution
}
\details{
\tabular{ll}{
Package: \tab KCsmart\cr
Type: \tab Package\cr
Version: \tab 2.3.2\cr
Date: \tab 2009-07-28\cr
License: \tab GPL\cr
}
Use the wrapper function 'calcSpm' to calculate the sample point matrix. Use 'findSigLevelTrad' to find a significance threshold using permutation based testing.
Use 'plot' to plot the sample point matrix or 'plotScaleSpace' to plot the significant regions over multiple scales (sigmas).
Use 'getSigSegments' to retrieve the significantly gained and lost regions using specific cutoffs.
To use the comparative version of KCsmart, use the 'calcSpmCollection', 'compareSpmCollection' and 'getSigRegionsCompKC' functions. See the documentation of those function for details on how to use these.
}
\author{
Jorma de Ronde, Christiaan Klijn

Maintainer: Jorma de Ronde <j.d.ronde@nki.nl>
}
\references{
Identification of cancer genes using a statistical framework for multiexperiment analysis of nondiscretized array CGH data.
Nucleic Acids Res. 2008 Feb;36(2):e13.
}

\keyword{ package }
\seealso{ 
\code{\link{calcSpm}},  \code{\link{findSigLevelTrad}}, \code{\link{findSigLevelFdr}}, \code{\link{plot}}, \code{\link{plotScaleSpace}}, \code{\link{getSigSegments}}
}
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)
spm4mb <- calcSpm(hsSampleData, hsMirrorLocs, sigma=4000000)

plot(spm1mb)
plot(spm1mb, chromosomes=c(1,5,6,'X'))

siglevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)
siglevel4mb <- findSigLevelTrad(hsSampleData, spm4mb, n=3)

plot(spm1mb, sig.level=siglevel1mb)

plotScaleSpace(list(spm1mb, spm4mb), list(siglevel1mb, siglevel4mb), type='g')

sigSegments1mb <- getSigSegments(spm1mb, siglevel1mb)


spmc1mb <- calcSpmCollection(hsSampleData, hsMirrorLocs, cl=c(rep(0,10),rep(1,10)))
spmcc1mb <- compareSpmCollection(spmc1mb, nperms=3)
spmcc1mbSigRegions <- getSigRegionsCompKC(spmcc1mb)

plot(spmcc1mb, sigRegions=spmcc1mbSigRegions)


}
