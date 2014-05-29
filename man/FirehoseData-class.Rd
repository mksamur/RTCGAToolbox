\name{FirehoseData-class}
\Rdversion{1.1}
\docType{class}
\alias{FirehoseData-class}

\title{Class \code{"FirehoseData"}}
\description{
This class is the main output of Firehose data client. All data will be downloaded as a FirehoseData class.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FirehoseData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Dataset}:}{Object of class \code{"character"}: Name of the TCGA dataset. }
    \item{\code{Clinical}:}{Object of class \code{"data.frame"}: Clinical dataset.}
    \item{\code{RNASeqGene}:}{Object of class \code{"matrix"}: Gene level RNASeq data. Data could be either raw counts or normalized.}
    \item{\code{RNASeq2GeneNorm}:}{Object of class \code{"matrix"}: Gene level RNASeq2 data. Firehose has 2 pre process pipeline for RNASeq data. More information can be accessible from Firehoseproject web site. Data could be either raw counts or normalized.}
    \item{\code{miRNASeqGene}:}{Object of class \code{"matrix"}: miRNA expression data from RNA seq platrom. Data could be either normalized or raw counts.}
    \item{\code{CNASNP}:}{Object of class \code{"data.frame"}: Segmented copy number altereration data from  SNP array.}
    \item{\code{CNVSNP}:}{Object of class \code{"data.frame"}: Segmented copy number variations data from  SNP array.}
    \item{\code{CNAseq}:}{Object of class \code{"data.frame"}: Segmented copy number alteration data from sequencing data.}
    \item{\code{CNACGH}:}{Object of class \code{"list"}: Lsit of segmented copy number alterations data from CGH array.}
    \item{\code{Methylation}:}{Object of class \code{"list"}: List of methylation array data. Each member of the list will be data from different platform.}
    \item{\code{mRNAArray}:}{Object of class \code{"list"}: List of mRNA expression data object. Each member of the list will be data from different platform. }
    \item{\code{miRNAArray}:}{Object of class \code{"list"}: List of miRNA expression data object. Each member of the list will be data from different platform. }
    \item{\code{RPPAArray}:}{Object of class \code{"list"}: List of RPPAArray data object. Each member of the list will be data from different platform. }
    \item{\code{Mutations}:}{Object of class \code{"data.frame"}: Mutation data related with samples.}
    \item{\code{GISTIC}:}{Object of class \code{"FirehoseGISTIC"}: Processed copy number estimations from GISTIC output of Firehose project.}
  }
}
\section{Methods}{
No methods defined with class "FirehoseData" in the signature.
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
\examples{
showClass("FirehoseData")
}
\keyword{classes}
