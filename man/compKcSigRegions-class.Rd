\name{compKcSigRegions-class}
\docType{class}
\alias{compKcSigRegions-class}
\alias{show,compKcSigRegions-method}

\title{KC smart comparative}
\description{	A matrix containing the results the significant regions for a given compKc object and FDR.}
\section{Objects from the Class}{
Objects can not be created by the user directly but rather through getSigRegionsCompKC.
}
\section{Slots}{
    \describe{
    \item{\code{regionTable}:}{The significant regions}
    \item{\code{method}:}{The method used to create the null distribution}
    \item{\code{cutoff}:}{The cutoff for the given false discovery rate which was used to determine the significant regions}
    \item{\code{fdr}:}{The false discovery rate used to determine the significant regions}

  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "compKcSigRegions")}: ... }
    \item{write.table}{\code{signature(object = "compKcSigRegions")}: ... }
	 }
}
\examples{
showClass("compKcSigRegions")
}
\keyword{classes}
