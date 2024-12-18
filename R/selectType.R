#' @include RTCGAToolbox-Class.R
NULL

#' Accessor function for the FirehoseData object
#'
#' An accessor function for the [RTCGAToolbox::FirehoseData-class]. An argument
#' will specify the data type to return See [RTCGAToolbox::FirehoseData-class]
#' for more details.
#'
#' @details
#' * clinical: Get the clinical data slot
#' * RNASeqGene: RNASeqGene
#' * RNASeq2GeneNorm: Normalized
#' * miRNASeqGene: micro RNA SeqGene
#' * CNASNP: Copy Number Alteration
#' * CNVSNP: Copy Number Variation
#' * CNASeq: Copy Number Alteration
#' * CNACGH: Copy Number Alteration
#' * Methylation: Methylation
#' * mRNAArray: Messenger RNA
#' * miRNAArray: micro RNA
#' * RPPAArray: Reverse Phase Protein Array
#' * Mutation: Mutations
#' * GISTIC: GISTIC v2 scores and probabilities
#'
#' @param object A `FirehoseData` class object
#' @param dataType A data type, see details.
#' @return The data type element of the `FirehoseData` object
setGeneric("selectType", function(object, dataType) standardGeneric("selectType"))


#' @describeIn FirehoseData Extract data type
#' @exportMethod selectType
#' @param dataType An available data type, see object show method
setMethod("selectType", "FirehoseData", function(object, dataType) {
    getElement(object, dataType)
})
