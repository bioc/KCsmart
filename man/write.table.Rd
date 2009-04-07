\name{write.table}
\docType{methods}
\alias{write.table}
\alias{write.table,sigSegments-method}
\alias{write.table,compKcSigRegions-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write summary of the significant regions to a table}
\description{
Write summary of the significant regions to a table
}
\usage{
write.table(x, file, append, quote, sep, eol, na, dec, row.names, col.names, qmethod)
\S4method{write.table}{sigSegments}(x, file="", append = FALSE, quote = 7, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names=c("#", "Status", "Chromosome", "Start", "End", "Average KC score", "Mode KC score", "Probes"), qmethod = c("escape", "double"))
\S4method{write.table}{compKcSigRegions}(x, file="", append = FALSE, quote = 7, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names=c("startrow", "endrow", "chromosome", "startposition", "endposition"), qmethod = c("escape", "double")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The sigSegments object to be summarized}
  \item{file}{ either a character string naming a file or a connection open
          for writing.  '""' indicates output to the console.}

  \item{append}{ logical. Only relevant if 'file' is a character string.  If
          'TRUE', the output is appended to the file.  If 'FALSE', any
          existing file of the name is destroyed.}

  \item{quote}{a logical value ('TRUE' or 'FALSE') or a numeric vector.  If
          'TRUE', any character or factor columns will be surrounded by
          double quotes.  If a numeric vector, its elements are taken
          as the indices of columns to quote.  In both cases, row and
          column names are quoted if they are written.  If 'FALSE',
          nothing is quoted.}

  \item{sep}{ the field separator string.  Values within each row of 'x'
          are separated by this string.}

  \item{eol}{ the character(s) to print at the end of each line (row).}

  \item{na}{ the string to use for missing values in the data.}

  \item{dec}{ the string to use for decimal points in numeric or complex
          columns: must be a single character.}

  \item{row.names}{ either a logical value indicating whether the row names of
          'x' are to be written along with 'x', or a character vector
          of row names to be written.}

  \item{col.names}{either a logical value indicating whether the column names
          of 'x' are to be written along with 'x', or a character
          vector of column names to be written.  See the section on
          'CSV files' for the meaning of 'col.names = NA'.}

   \item{qmethod}{ a character string specifying how to deal with embedded
          double quote characters when quoting strings.  Must be one of
          '"escape"' (default), in which case the quote character is
          escaped in C style by a backslash, or '"double"', in which
          case it is doubled.  You can specify just the initial letter.}


}
\details{
Writes a summary of the sigSegments object to file. The resulting table contains 7 columns. The interpretation of the columns is as follows:
    \item{Status}{Either 'L' for loss or 'G' for gain}
    \item{Chromosome}{The chromosome on which this segment is located}
    \item{Start}{The start position (in base pairs) of the segment on the chromosome}
    \item{End}{The end position of the segment on the chromosome}
    \item{Average KC score}{The average KCsmart score over all base pairs in this segment}
    \item{Mode KC score}{The highest (for gains) or lowest (for losses) KCsmart score over all base pairs in this segment}
    \item{Probes}{All probes from the original data that fall into this segment}
}
\value{

}
\author{ Jorma de Ronde }
\note{ 
}
\seealso{  \code{\link{calcSpm}}, \code{\link{getSigSegments}} }
\examples{
data(hsSampleData)
data(hsMirrorLocs)

spm1mb <- calcSpm(hsSampleData, hsMirrorLocs)

siglevel1mb <- findSigLevelTrad(hsSampleData, spm1mb, n=3)

sigSegments1mb <- getSigSegments(spm1mb, siglevel1mb)
write.table(sigSegments1mb, file=file.path(tempdir(),'sigSegments1mb.txt'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}
