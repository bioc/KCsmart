\name{samplePointMatrix-class}
\docType{class}
\alias{samplePointMatrix-class}
\alias{show,samplePointMatrix-method}

\title{Sample point matrix}
\description{	A sample point matrix resulting from a call to calcSpm}
\section{Objects from the Class}{
Objects can not be created by the user directly but rather through calcSpm.
}
\section{Slots}{
    \describe{
    \item{\code{totalLength}:}{ The total length of the sample point matrix, measures in sample points}
    \item{\code{maxy}:}{The maximum KC score attained over the sample point matrix }
    \item{\code{miny}:}{The minimum KC score attained over the sample point matrix }
    \item{\code{sampleDensity}:}{The sample density used to calculate the sample point matrix. ie the distance between two points in the sample point matrix, measured in base pairs. }
    \item{\code{sigma}:}{The sigma used for the kernel to calculate the sample point matrix.}
    \item{\code{mirrorLocs}:}{The mirror locations list used to calculate the sample point matrix}
    \item{\code{probeAnnotation}:}{The original probe annotation from the input data.}
    \item{\code{data}:}{The sample point matrix data points in the form of a list where each list element represents a chromosome.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "samplePointMatrix")}: ... }
    \item{show}{\code{signature(object = "samplePointMatrix")}: ... }
	 }
}
\examples{
showClass("samplePointMatrix")
}
\keyword{classes}
