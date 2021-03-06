% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFirehoseData.R
\name{getFirehoseData}
\alias{getFirehoseData}
\title{Get data from Firehose portal.}
\usage{
getFirehoseData(
  dataset,
  runDate = "20160128",
  gistic2Date = "20160128",
  RNASeqGene = FALSE,
  RNASeq2Gene = FALSE,
  clinical = TRUE,
  miRNASeqGene = FALSE,
  RNASeq2GeneNorm = FALSE,
  CNASNP = FALSE,
  CNVSNP = FALSE,
  CNASeq = FALSE,
  CNACGH = FALSE,
  Methylation = FALSE,
  Mutation = FALSE,
  mRNAArray = FALSE,
  miRNAArray = FALSE,
  RPPAArray = FALSE,
  GISTIC = FALSE,
  RNAseqNorm = "raw_count",
  RNAseq2Norm = "normalized_count",
  forceDownload = FALSE,
  destdir = tempdir(),
  fileSizeLimit = 500,
  getUUIDs = FALSE,
  ...
)
}
\arguments{
\item{dataset}{A cohort name. TCGA cancer code obtained via \code{\link{getFirehoseDatasets}}}

\item{runDate}{Standard data run dates. Date list can be accessible via \code{\link{getFirehoseRunningDates}}}

\item{gistic2Date}{Analysis run date for GISTIC obtained via \code{\link{getFirehoseAnalyzeDates}}}

\item{RNASeqGene}{Logical (default FALSE) RNAseq TPM data.}

\item{RNASeq2Gene}{Logical (default FALSE) RNAseq v2 (RSEM processed) data; see \code{RNAseqNorm} argument.}

\item{clinical}{Logical (default TRUE) clinical data.}

\item{miRNASeqGene}{Logical (default FALSE) smallRNAseq data.}

\item{RNASeq2GeneNorm}{Logical (default FALSE) RNAseq v2 (RSEM processed) data.}

\item{CNASNP}{Logical (default FALSE) somatic copy number alterations data from SNP array.}

\item{CNVSNP}{Logical (default FALSE) germline copy number variants data from SNP array.}

\item{CNASeq}{Logical (default FALSE) somatic copy number alterations data from sequencing.}

\item{CNACGH}{Logical (default FALSE) somatic copy number alterations data from CGH.}

\item{Methylation}{Logical (default FALSE) methylation data.}

\item{Mutation}{Logical (default FALSE) mutation data from sequencing.}

\item{mRNAArray}{Logical (default FALSE) mRNA expression data from microarray.}

\item{miRNAArray}{Logical (default FALSE) miRNA expression data from microarray.}

\item{RPPAArray}{Logical (default FALSE) RPPA data}

\item{GISTIC}{logical (default FALSE) processed copy number data}

\item{RNAseqNorm}{RNAseq data normalization method. (Default raw_count)}

\item{RNAseq2Norm}{RNAseq v2 data normalization method. (Default normalized_count, raw_count, scaled_estimate)}

\item{forceDownload}{A logic (Default FALSE) key to force download RTCGAToolbox every time. By default if you download files into your working directory once than RTCGAToolbox using local files next time.}

\item{destdir}{Directory in which to store the resulting downloaded file.
Defaults to a temporary directory given by \code{tempdir()}.}

\item{fileSizeLimit}{Files that are larger than set value (megabyte) won't be downloaded (Default: 500)}

\item{getUUIDs}{Logical key to get UUIDs from barcode (Default: FALSE)}

\item{...}{Additional arguments to pass down.}
}
\value{
A \code{FirehoseData} data object that stores data for selected data types.
}
\description{
\code{getFirehoseData} returns \code{FirehoseData} object that stores TCGA data.
}
\details{
This is a main client function to download data from Firehose TCGA portal.
}
\examples{
# Sample Dataset
data(accmini)
accmini
\dontrun{
BRCAdata <- getFirehoseData(dataset="BRCA",
runDate="20140416",gistic2Date="20140115",
RNASeqGene=TRUE,clinical=TRUE,mRNAArray=TRUE,Mutation=TRUE)
}
}
\seealso{
\link{getLinks}
}
