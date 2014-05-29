\name{getDiffExpressedGenes}
\alias{getDiffExpressedGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get differentially expressed genes between tumor and normal samples.
%%  ~~function to do ... ~~
}
\description{
This function inherit differential gene expression analysis from "limma" package. The function creates 2 groups (normal and tumor) based on TCGA barcodes and calculate differential gene expression between groups. Finally it returns "DGEResult" object list.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getDiffExpressedGenes(dataObject, DrawPlots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataObject}{: "FirehoseData" object that includes raw gene count from RNAseq platform or normalized expression values from array platfroms.
%%     ~~Describe \code{dataObject} here~~
}
  \item{DrawPlots}{: if set as "TRUE", function will draw heatmap, volcano plot and model plot from "voom" function.
%%     ~~Describe \code{DrawPlots} here~~
}
}
\value{
This function returns a list where each item will be "DGEResult" object for one comparison.
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
data(RTCGASample)
dgegenes = getDiffExpressedGenes(a2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
