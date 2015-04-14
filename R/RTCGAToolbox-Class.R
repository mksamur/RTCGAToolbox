#' An S4 class to store data from CGA platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data frame that stores the CGH data.
setClass("FirehoseCGHArray", representation(Filename = "character", DataMatrix = "data.frame"))

#' An S4 class to store data from methylation platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data frame that stores the methylation data.
setClass("FirehoseMethylationArray", representation(Filename = "character", DataMatrix = "data.frame"))

#' An S4 class to store data from array (mRNA, miRNA etc.) platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data matrix that stores the expression data.
setClass("FirehosemRNAArray", representation(Filename = "character", DataMatrix = "matrix"))

#' An S4 class to store processed copy number data. (Data processed by using GISTIC2 algorithm)
#'
#' @slot Dataset Cohort name
#' @slot AllByGene A data frame that stores continuous copy number
#' @slot ThresholedByGene A data frame for discrete copy number data
setClass("FirehoseGISTIC", representation(Dataset = "character", AllByGene = "data.frame",ThresholedByGene="data.frame"))

#' An S4 class to store main data object from clinent function.
#'
#' @slot Dataset A cohort name
#' @slot Clinical Clinical data frame
#' @slot RNASeqGene Gene level expression data matrix from RNAseq
#' @slot RNASeq2GeneNorm Gene level expression data matrix from RNAseq (RSEM)
#' @slot miRNASeqGene miRNA expression data from matrix smallRNAseq
#' @slot CNASNP A data frame to store somatic copy number alterations from SNP array platform
#' @slot CNVSNP A data frame to store germline copy number variants from SNP array platform
#' @slot CNAseq A data frame to store somatic copy number alterations from sequencing platform
#' @slot CNACGH A list that stores \code{FirehoseCGHArray} object for somatic copy number alterations from CGH platform
#' @slot Methylation A list that stores \code{FirehoseMethylationArray} object for methylation data
#' @slot mRNAArray A list that stores \code{FirehosemRNAArray} object for gene expression data from microarray
#' @slot miRNAArray A list that stores \code{FirehosemRNAArray} object for miRNA expression data from microarray
#' @slot RPPAArray A list that stores \code{FirehosemRNAArray} object for RPPA data
#' @slot Mutations A data frame for mutation infromation from sequencing data
#' @slot GISTIC A \code{FirehoseGISTIC} object to store processed copy number data
#' @slot BarcodeUUID A data frame that stores the Barcodes, UUIDs and Short sample identifiers
setClass("FirehoseData", representation(Dataset = "character", Clinical = "data.frame", RNASeqGene = "matrix",
                                        RNASeq2GeneNorm="matrix",miRNASeqGene="matrix",CNASNP="data.frame",
                                        CNVSNP="data.frame",CNAseq="data.frame",CNACGH="list",Methylation="list",
                                        mRNAArray="list",miRNAArray="list",RPPAArray="list",Mutations="data.frame",
                                        GISTIC="FirehoseGISTIC",BarcodeUUID="data.frame"))
setMethod("show", "FirehoseData",function(object){
  message(paste0(object@Dataset," FirehoseData object"))
  message("Available slots:")
  if(dim(object@Clinical)[1] > 0 & dim(object@Clinical)[2] > 0){message("@Clinical: A data frame, dim: ",paste(dim(object@Clinical),collapse = "\t"))}
  if(dim(object@RNASeqGene)[1] > 0 & dim(object@RNASeqGene)[2] > 0){message("@RNASeqGene: A matrix withraw read counts or normalized data, dim: ",paste(dim(object@RNASeqGene),collapse = "\t"))}
  if(dim(object@RNASeq2GeneNorm)[1] > 0 & dim(object@RNASeq2GeneNorm)[2] > 0){message("@RNASeq2GeneNorm: A matrix withraw read counts or normalized data, dim: ",paste(dim(object@RNASeq2GeneNorm),collapse = "\t"))}
  if(dim(object@miRNASeqGene)[1] > 0 & dim(object@miRNASeqGene)[2] > 0){message("@miRNASeqGene: A matrix, dim: ",paste(dim(object@miRNASeqGene),collapse = "\t"))}
  if(dim(object@CNASNP)[1] > 0 & dim(object@CNASNP)[2] > 0){message("@CNASNP: A data.frame, dim: ",paste(dim(object@CNASNP),collapse = "\t"))}
  if(dim(object@CNVSNP)[1] > 0 & dim(object@CNVSNP)[2] > 0){message("@CNVSNP: A data.frame, dim: ",paste(dim(object@CNVSNP),collapse = "\t"))}
  if(dim(object@CNAseq)[1] > 0 & dim(object@CNAseq)[2] > 0){message("@CNAseq: A data.frame, dim: ",paste(dim(object@CNAseq),collapse = "\t"))}
  if(length(object@CNACGH) > 0 ){message("@CNACGH: A list contains FirehoseCGHArray object(s), length: ",length(object@CNACGH))}
  if(length(object@Methylation) > 0 ){message("@Methylation: A list contains FirehoseMethylationArray object(s), length: ",length(object@Methylation))}
  if(length(object@mRNAArray) > 0 ){message("@mRNAArray: A list contains FirehosemRNAArray object(s), length: ",length(object@mRNAArray))}
  if(length(object@miRNAArray) > 0 ){message("@miRNAArray: A list contains FirehosemRNAArray object(s), length: ",length(object@miRNAArray))}
  if(length(object@RPPAArray) > 0 ){message("@RPPAArray: A list contains FirehosemRNAArray object(s), length: ",length(object@RPPAArray))}
  if(length(object@GISTIC@Dataset) > 0){message("@GISTIC: A FirehoseGISTIC object to store copy number data")}
  if(dim(object@Mutations)[1] > 0 & dim(object@Mutations)[2] > 0){message("@Mutations: A data.frame, dim: ",paste(dim(object@Mutations),collapse = "\t"))}
}
)

#' An S4 class to store differential gene expression results
#'
#' @slot Dataset Dataset name
#' @slot Toptable Results data frame
setClass("DGEResult", representation(Dataset = "character", Toptable = "data.frame"))

#' An S4 class to store correlations between gene expression level and copy number data
#'
#' @slot Dataset A cohort name
#' @slot Correlations Results data frame
setClass("CorResult", representation(Dataset = "character", Correlations = "data.frame"))