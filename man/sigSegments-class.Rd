\name{sigSegments-class}
\docType{class}
\alias{sigSegments-class}
\alias{show,sigSegments-method}

\title{Significant segments }
\description{Lists the significant segments found in a given sample point matrix using a given significance level}
\section{Objects from the Class}{
Objects can not be created by the user directly but rather through getSigSegments.
}
\section{Slots}{
	 \describe{
    \item{\code{gains}:}{Gained segments}
    \item{\code{losses}:}{Lost segments}
    \item{\code{sigma}:}{The sigma used for the kernel to calculate the sample point matrix.}
    \item{\code{sigLevels}:}{The significance levels at which significant segments are calculated}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "sigSegments")}: ... }
    \item{write.table}{\code{signature(x = "sigSegments")}: ... }
	 }
}
\examples{
showClass("sigSegments")
}
\keyword{classes}
