#' @include RTCGAToolbox-Class.R
NULL

#' @name selectType
#' @title Accessor function for the FirehoseData object
#' @description An accessor function for the \linkS4class{FirehoseData}
#' class. An argument will specify the data type to return
#' See \link{FirehoseData-class} for more details.
#'
#' @details
#' \itemize{
#'     \item{clinical} - Get the clinical data slot
#'     \item{RNASeqGene} - RNASeqGene
#'     \item{RNASeq2GeneNorm} - Normalized
#'     \item{miRNASeqGene} - micro RNA SeqGene
#'     \item{CNASNP} - Copy Number Alteration
#'     \item{CNVSNP} - Copy Number Variation
#'     \item{CNASeq} - Copy Number Alteration
#'     \item{CNACGH} - Copy Number Alteration
#'     \item{Methylation} - Methylation
#'     \item{mRNAArray} - Messenger RNA
#'     \item{miRNAArray} - micro RNA
#'     \item{RPPAArray} - Reverse Phase Protein Array
#'     \item{Mutation} - Mutations
#'     \item{GISTIC} - GISTIC v2 scores and probabilities
#' }
#'
#' @param object A \code{FirehoseData} class object
#' @param dataType A data type, see details.
#' @return The data type element of the \code{FirehoseData} object
setGeneric("selectType", function(object, dataType) standardGeneric("selectType"))


#' @describeIn FirehoseData Extract data type
#' @aliases NULL
#' @exportMethod selectType
setMethod("selectType", "FirehoseData", function(object, dataType) {
    getElement(object, dataType)
})
