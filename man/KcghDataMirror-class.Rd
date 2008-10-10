\name{KcghDataMirror-class}
\docType{class}
\alias{KcghDataMirror-class}
\alias{initialize,KcghDataMirror-method}

\title{Class "KcghDataMirror"}
\description{Internal class}

\section{Slots}{
    \describe{
    \item{\code{mirrorLocs}:}{Holds mirrorLocs object }
    \item{\code{probeAnnotation}:}{Object of class \code{"probeAnnotation"}}
    \item{\code{pos}:}{Holds aCGH data for losses}
    \item{\code{neg}:}{Holds aCGH data for gains}
    \item{\code{nrSamples}:}{The number of samples in this analysis}
  }
}
\section{Extends}{
Class \code{"\linkS4class{KcghDataSum}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "KcghDataMirror")}: For internal use only }
	 }
}
\author{Jorma de Ronde }
\note{For internal use only }

\keyword{classes}
