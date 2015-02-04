\name{getFirehoseData}
\alias{getFirehoseData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Firehose project data client
%%  ~~function to do ... ~~
}
\description{
This function will download data from Firehose project portal for further steps. Downloaded data will be stored as a "FirehoseData" object and text files will be stored in working directory by default.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getFirehoseData(dataset, runDate = NULL, gistic2_Date = NULL, 
RNAseq_Gene = FALSE, Clinic = TRUE, miRNASeq_Gene = FALSE, 
RNAseq2_Gene_Norm = FALSE, CNA_SNP = FALSE, CNV_SNP = FALSE, 
CNA_Seq = FALSE,CNA_CGH = FALSE, Methylation = FALSE, 
Mutation = FALSE, mRNA_Array = FALSE,miRNA_Array = FALSE, 
RPPA = FALSE,RNAseqNorm = "raw_counts", RNAseq2Norm = "normalized_count", todir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{: Dataset TCGA alias. A valid list can be accessible with "getFirehoseDatasets" funtion. 
%%     ~~Describe \code{dataset} here~~
}
  \item{runDate}{: Firehose data run date. A valid dates can be accessible with "getFirehoseRunningDates" function.
%%     ~~Describe \code{runDate} here~~
}
  \item{gistic2_Date}{: Firehose data analyze run date. A valid dates can be accessible with "getFirehoseAnalyzeDates" function.
%%     ~~Describe \code{gistic2_Date} here~~
}
  \item{RNAseq_Gene}{: If "TRUE" client will donwload gene lebel RNASeq data.
%%     ~~Describe \code{RNAseq_Gene} here~~
}
  \item{Clinic}{: If "TRUE" client will download clinical data.
%%     ~~Describe \code{Clinic} here~~
}
  \item{miRNASeq_Gene}{: If "TRUE" client will download miRNA gene level data.
%%     ~~Describe \code{miRNASeq_Gene} here~~
}
  \item{RNAseq2_Gene_Norm}{: If "TRUE" client will download gene level RNASeq data from second pipeline.
%%     ~~Describe \code{RNAseq2_Gene_Norm} here~~
}
  \item{CNA_SNP}{: If "TRUE" client will download client will download copy number alterations data processed by SNP arrays.
%%     ~~Describe \code{CNA_SNP} here~~
}
  \item{CNV_SNP}{: If "TRUE" client will download client will download copy number variation data processed by SNP arrays.
%%     ~~Describe \code{CNV_SNP} here~~
}
  \item{CNA_Seq}{: If "TRUE" client will download client will download copy number alterations data processed by DNA sequencing.
%%     ~~Describe \code{CNA_Seq} here~~
}
  \item{CNA_CGH}{: If "TRUE" client will download client will download copy number alterations data processed by CGH arrays.
%%     ~~Describe \code{CNA_CGH} here~~
}
  \item{Methylation}{: If "TRUE" client will download methylation array data.
%%     ~~Describe \code{Methylation} here~~
}
  \item{Mutation}{: If "TRUE" client will download mutation information.
%%     ~~Describe \code{Mutation} here~~
}
  \item{mRNA_Array}{: If "TRUE" client will download gene level expression profile processed by microarrays.
%%     ~~Describe \code{mRNA_Array} here~~
}
  \item{miRNA_Array}{: If "TRUE" client will download miRNA gene level expression profile processed by microarrays.
%%     ~~Describe \code{miRNA_Array} here~~
}
  \item{RPPA}{: If "TRUE" client will download reverse phase protein array data.
%%     ~~Describe \code{RPPA} here~~
}
  \item{RNAseqNorm}{: RNAseq data type normalization. Default: "raw_counts". Users can set "RPKM" depends on further analysis.
%%     ~~Describe \code{RNAseqNorm} here~~
}
  \item{RNAseq2Norm}{: RNAseq2 data type normalization. Default: "normalized_count". This method uses different algorithms to determine the gene expression levels. Users can set "RPKM" depends on further analysis.
%%     ~~Describe \code{RNAseq2Norm} here~~
}
  \item{todir}{: The directory to download files to. It will be created if necessary. The default is to use the current working directory.  
%%     ~~Describe \code{todir} here~~
}
}
\value{
This function returns a "FirehoseData" object.
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
\dontrun{a2 <- getFirehoseData(dataset="BRCA",
                  runDate="20140416",gistic2_Date="20140115",
                  RNAseq_Gene=TRUE,Clinic=TRUE,mRNA_Array=TRUE,Mutation=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
