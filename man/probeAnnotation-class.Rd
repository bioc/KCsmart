\name{probeAnnotation-class}
\docType{class}
\alias{probeAnnotation-class}
\alias{[,probeAnnotation-method}
\alias{initialize,probeAnnotation-method}

\title{Class "probeAnnotation"}
\description{Holds the probe annotation}
\section{Objects from the Class}{
 Instances of this class are not meant to be created by the user
}

\section{Slots}{
	 \describe{
    \item{\code{chromosome}:}{Chromosome on which the probe is located}
    \item{\code{maploc}:}{Location of the probe on the chromosome}
    \item{\code{name}:}{Probe name}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "probeAnnotation")}: Access information about a probe }
    \item{initialize}{\code{signature(.Object = "probeAnnotation")}: Internal use only }
	 }
}

\author{Jorma de Ronde }

\keyword{classes}
