\name{spmCollection-class}
\docType{class}
\alias{spmCollection-class}
\alias{show,spmCollection-method}

\title{Sample point matrix collection}
\description{	A sample point matrix collection resulting from a call to calcSpmCollection}
\section{Objects from the Class}{
Objects can not be created by the user directly but rather through calcSpmCollection.
}
\section{Slots}{
    \describe{
    \item{\code{annotation}:}{The annotation (containing the chromosome and position on the chromosome) for the sample points in the 'data' slot}
    \item{\code{data}:}{A matrix with the kernel smoothed estimates of all samples}
    \item{\code{cl}:}{A class vector indicating which samples belong to which class}
    \item{\code{sampleDensity}:}{The sample density used to calculate the sample point matrix. ie the distance between two points in the sample point matrix, measured in base pairs. }
    \item{\code{sigma}:}{The sigma used for the kernel to calculate the sample point matrix.}
    \item{\code{mirrorLocs}:}{The mirror locations list used to calculate the sample point matrix}

  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "spmCollection")}: ... }
	 }
}
\examples{
showClass("spmCollection")
}
\keyword{classes}
