\name{CorResult-class}
\Rdversion{1.1}
\docType{class}
\alias{CorResult-class}

\title{Class \code{"CorResult"}}
\description{
This object stores correlation analysis results between copy number and gene expression data.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CorResult", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Dataset}:}{Object of class \code{"character"}: Stores name of the dataset. }
    \item{\code{Correlations}:}{Object of class \code{"data.frame"}: Stores correlation analysis results. First colum will be gene symbol, second column will be correlation value and last column will be BH adjusted p value.}
  }
}
\section{Methods}{
No methods defined with class "CorResult" in the signature.
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
\examples{
showClass("CorResult")
}
\keyword{classes}
