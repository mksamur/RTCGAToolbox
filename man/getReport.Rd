\name{getReport}
\alias{getReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make report plot.
%%  ~~function to do ... ~~
}
\description{
This function provides a circle plot from "FirehoseData" data object. Current function reports differentially expressed genes, copy number values and mutations. Function can take 2 differentially gene expressin analysis results. Users may provide results from 2 different platform. Output will be saved as apdf file.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getReport(dataObject, DGEResult1 = NULL, DGEResult2 = NULL, geneLocations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataObject}{: "FirehoseData" object that includes copy number data and mutation data.
%%     ~~Describe \code{dataObject} here~~
}
  \item{DGEResult1}{: Differentially gene expression results. Object type must be "DGEResult". Report plot will use fold change value to highlight up and down regulated genes.
%%     ~~Describe \code{DGEResult1} here~~
}
  \item{DGEResult2}{: Differentially gene expression results. Object type must be "DGEResult". Report plot will use fold change value to highlight up and down regulated genes.
%%     ~~Describe \code{DGEResult2} here~~
}
  \item{geneLocations}{: Gene location data frame. RTCGAToolbox provides location for hg19 known genes from UCSC browser. Plese check "RTCGASample" data.
%%     ~~Describe \code{geneLocations} here~~
}
}
\value{
This function creates a pdf file on the working directory that shows the differentially expressed genes, copy number and mutation information.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
\examples{
\dontrun{
data(RTCGASample)
data(hg19.ucsc.gene.locations)
t1=getDiffExpressedGenes(a2)
getReport(dataObject=a2,DGEResult1=t1[[1]],geneLocations=hg19.ucsc.gene.locations)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
