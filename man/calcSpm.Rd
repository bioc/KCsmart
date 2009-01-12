\name{calcSpm}
\alias{calcSpm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ KCsmart wrapper }
\description{
  Wrapper function that calculates the sample point matrix from the aCGH data
}
\usage{
calcSpm(data, mirrorLocs, sigma = 1e+06, sampleDensity = 50000, maxmem = 1000, verbose=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The aCGH data. Can either be in DNAcopy format or as a data.frame described in the details section }
  \item{mirrorLocs}{ List containing the chromosome start, centromere and end positions }
  \item{sigma}{ The kernel width }
  \item{sampleDensity}{ The sample point matrix resolution }
  \item{maxmem}{ This parameter controls memory usage, set to lower value to lower memory consumption }
  \item{verbose}{ If set to false, no progress information is displayed }
}
\details{
  'data' can be in cghRaw (CGHbase), DNAcopy or in data.frame format. When using the latter, the data.frame must have the following two columns: 'chrom' stating the chromosome
  the probe is located on, 'maploc' describing the position on the chromosome of the probe. The remainder of the data.frame will be interpreted as sample data points.
  The row names of that data will be used as probe names (when available).
Important note: the data can not contain any missing values. If your data includes missing values you will need to preprocess (for example impute) it using other software solutions.
  
  The mirror locations for Homo Sapiens and Mus Musculus are provided in the package. These can be loaded using data(hsMirrorLocs) and data(mmMirrorLocs) 
  respectively. The 'mirrorLocs' object is a list with vectors containing the start, centromere (optional) and end of each chromosome as the list elements. Additionally it should 
  contain an attribute 'chromNames' listing the chromosome names of each respective list element.
  
  'sigma' defines the kernel width of the kernel used to convolute the data.
  
  'sampleDensity' defines the resolution of the sample point matrix to be calculated. A sampleDensity of 50000 would correspond to a sample point every 50k base pairs.
  
}
\value{
 Returns a sample point matrix object. The object has several slots of which the 'data' slot contains a list where each list item represents a chromosome. Each list item in turn contains the sample point matrix
  for the gains and the losses separately and an attribute specifying the corresponding chromosome.
  The sample point matrix contains the following additional slots:
  totalLength: Total length of the sample point matrix
  maxy and miny: Maximal and minimal score attained
  
  The other slots just represent the parameters used to calculate the sample point matrix.
  
  Use 'plot' to plot the sample point matrix and 'findSigLevelTrad' to find a significance threshold.
  'plotScaleSpace' can be used to plot the significant regions of multiple sample point matrices (using different sigmas).
}
\references{ }
\author{ Jorma de Ronde }
\note{ 
}
\seealso{ 
\code{\link{plot}},  \code{\link{findSigLevelTrad}}, \code{\link{plotScaleSpace}}
}
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)
spm4mb <- calcSpm(hsSampleData, hsMirrorLocs, sigma=4000000)

plot(spm1mb)
plot(spm1mb, chromosomes=c(1,5,6,'X'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
