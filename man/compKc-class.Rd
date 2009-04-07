\name{compKc-class}
\docType{class}
\alias{compKc-class}
\alias{show,compKc-method}
\alias{initialize,compKc-method}
\alias{plot,compKc-method}

\title{KC smart comparative}
\description{	A matrix containing the results from a call to compareSpmCollection}
\section{Objects from the Class}{
Objects can not be created by the user directly but rather through compareSpmCollection.
}
\section{Slots}{
    \describe{
    \item{\code{spmCollection}:}{The original spmCollection used to compare the samples}
    \item{\code{method}:}{The method used to create the null distribution}
    \item{\code{siggenesResult}:}{In case of the siggenes method being used, a siggenes object containing the fdr-cutoff table}
    \item{\code{snrResult}:}{In case of the signal-noise/permutation based method being used, the signal-to-noise data and a matrix with the (class based) permutations}

  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "compKc")}: Internal use only }
  	\item{plot}{\code{signature(x = "compKc")}: ... }
    \item{show}{\code{signature(object = "compKc")}: ... }
	 }
}
\examples{
showClass("compKc")
}
\keyword{classes}
