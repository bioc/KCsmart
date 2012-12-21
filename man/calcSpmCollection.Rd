\name{calcSpmCollection}
\alias{calcSpmCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ KCsmart Comparative wrapper }
\description{
  Wrapper function that calculates the sample point matrix collection from the aCGH data. The sample point matrix collection is used in the comparative version of KCsmart.
}
\usage{
calcSpmCollection(data, mirrorLocs, cl=NULL, data2=NULL, sigma=1000000, sampleDensity=50000, maxmem=1000, verbose=F, doChecks=T, old=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The aCGH data. Can either be in DNAcopy format or as a data.frame described in the details section }
  \item{mirrorLocs}{ List containing the chromosome start, centromere and end positions }
  \item{cl}{ A class vector indicating which samples belong to which class}
  \item{data2}{ Instead of a class vector a second data set can be provided which will be combined with the first data set into one sample point matrix collection}
  \item{sigma}{ The kernel width }
  \item{sampleDensity}{ The sample point matrix resolution }
  \item{maxmem}{ This parameter controls memory usage, set to lower value to lower memory consumption }
  \item{verbose}{ If set to false, no progress information is displayed }
  \item{doChecks}{ If set to false, the data will not be checked for consistency }
  \item{old}{ If set to true the old implementation of KCsmart will be used to calculate the spm }
}
\details{
  The input can either consist of a single data set and a class vector or two separate datasets. In the latter case a class vector will be created assigning each data set to its own class.
  'data' can be in cghRaw (CGHbase), DNAcopy or in data.frame format. When using the latter, the data.frame must have the following two columns: 'chrom' stating the chromosome
  the probe is located on, 'maploc' describing the position on the chromosome of the probe. The remainder of the data.frame will be interpreted as sample data points.
  The row names of that data will be used as probe names (when available).
Important note: the data can not contain any missing values. If your data includes missing values you will need to preprocess (for example impute) it using other software solutions.
  
  The mirror locations for Homo Sapiens and Mus Musculus are provided in the package. These can be loaded using data(hsMirrorLocs) and data(mmMirrorLocs) 
  respectively. The 'mirrorLocs' object is a list with vectors containing the start, centromere (optional) and end of each chromosome as the list elements. Additionally it should 
  contain an attribute 'chromNames' listing the chromosome names of each respective list element.
  
  'sigma' defines the kernel width of the kernel used to convolute the data.
  
  'sampleDensity' defines the resolution of the sample point matrix to be calculated. A sampleDensity of 50000 would correspond to a sample point every 50k base pairs.
  
  'old' can be used if you want to reproduce data that was generated with old (pre 2.9.0) versions of KCsmart, for any new analyses we recommend this flag to be set to false
  
}
\value{
 Returns a sample point matrix collection object. The object has several slots of which the 'data' slot contains a matrix with the kernel smoothed estimates of all samples.
  The sample point matrix collection contains the following additional slots:
  cl: A class vector indicating which samples belong to which class. 
  annotation: The annotation (containing the chromosome and position on the chromosome) for the sample points in the 'data' slot
  
  The other slots just represent the parameters used to calculate the sample point matrix collection.
  
  Use 'compareSpmCollection' to get a 'compKc' object for which the significant regions can be calculated using 'getSigRegionsCompKC'.
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
