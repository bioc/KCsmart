\name{KcghDataSum-class}
\docType{class}
\alias{KcghDataSum-class}
\alias{initialize,KcghDataSum-method}
\alias{sort,KcghDataSum,missing-method}

\title{Class "KcghDataSum" }
\description{Internal class}


\section{Slots}{
    \describe{
    \item{\code{probeAnnotation}:}{Object of class \code{"probeAnnotation"}}
    \item{\code{pos}:}{Holds aCGH data for losses}
    \item{\code{neg}:}{Holds aCGH data for gains}
    \item{\code{nrSamples}:}{The number of samples in this analysis}
  }
}

\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "KcghDataSum")}: For internal use only}
    \item{sort}{\code{signature(x = "KcghDataSum", decreasing = "missing")}: For internal use only}
	 }
}

\author{Jorma de Ronde }
\note{For internal use only }

\keyword{classes}
