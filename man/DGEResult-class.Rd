\name{DGEResult-class}
\Rdversion{1.1}
\docType{class}
\alias{DGEResult-class}

\title{Class \code{"DGEResult"}}
\description{
This class stores differential gene expression analysis results.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DGEResult", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Dataset}:}{Object of class \code{"character"}: Name of the dataset eg:RNASeq or the chip platfrom}
    \item{\code{Toptable}:}{Object of class \code{"data.frame"}: Top table output from limma package.}
  }
}
\section{Methods}{
No methods defined with class "DGEResult" in the signature.
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
\examples{
showClass("DGEResult")
}
\keyword{classes}
