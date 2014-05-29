\name{FirehoseGISTIC-class}
\Rdversion{1.1}
\docType{class}
\alias{FirehoseGISTIC-class}

\title{Class \code{"FirehoseGISTIC"}}
\description{
This class stores the copy number estimation from GISTIC2 algorithm.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FirehoseGISTIC", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Dataset}:}{Object of class \code{"character"}: Dataset name}
    \item{\code{AllByGene}:}{Object of class \code{"data.frame"}: Gene level copy number estimation. 0 accepted as normal copy number. Negative values indicate deletions, positive values indicate amplifications.}
    \item{\code{ThresholedByGene}:}{Object of class \code{"data.frame"}: Thresholed copy number estimations.-2: possibly homozygous deletion, -1: possibly heterozygous deletion, 0: diploid, 1: gain, 2: high amplification  }
  }
}
\section{Methods}{
No methods defined with class "FirehoseGISTIC" in the signature.
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
\examples{
showClass("FirehoseGISTIC")
}
\keyword{classes}
