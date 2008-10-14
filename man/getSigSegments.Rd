\name{getSigSegments}
\alias{getSigSegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve the significantly gained and lost regions including the corresponding, original probes }
\description{
  Retrieve the significantly gained and lost regions including the corresponding, original probes. A significance level must be selected by the user.
}
\usage{
getSigSegments(spm, sigLevels, chromosomes=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{spm}{ The sample point matrix to be plotted }
  \item{sigLevels}{ The significance thresholds to be used}
  \item{chromosomes}{ Takes a vector of chromosomes to be plotted. Defaults to all chromosomes. }
}
\details{
  'sigLevels' should contain the significance thresholds in a list with the positive (gains) threshold in the 'pos' element and the negative (losses) threshold in the 'neg' element.
  This is the format as returned by 'findSigLevelTrad' and 'findSigLevelFdr'.
}
\value{
  Returns a sigSegments object containing the chromosome, start position, end position, average KC score and the mode of the KC score in that region of all segments passing the thresholds as set in 'sigLevels'. 
  Additionally, returns the IDs and indices of the probes and the positions in the sample point matrix within the significant regions.
  The results are stored in two separate slots: 'gains' for gains and 'losses' for losses.
  Use 'write.table' to save the results to file.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Jorma de Ronde }
\note{
}
\seealso{  \code{\link{findSigLevelTrad}},  \code{\link{findSigLevelTrad}}, \code{\link{write.table}}}
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)

siglevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)

sigSegments1mb <- getSigSegments(spm1mb, siglevel1mb)
\dontrun{write.table(sigSegments1mb, file=paste(tempdir(),'/sigSegments1mb.txt',sep=""))}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
