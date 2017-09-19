#' @include RTCGAToolbox-Class.R
NULL

#' @name FirehoseData-extractors
#' @title Accessor methods for the FirehoseData object
#' @description A list of accessor functions for the \linkS4class{FirehoseData}
#' class. Each exposing a unique element of the \linkS4class{FirehoseData}
#' object. See \link{FirehoseData-class} for more details.
#' 
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
#' @param object A \code{FirehoseData} class object
NULL


#' @rdname FirehoseData-extractors 
#' @export
setGeneric("clinical", function(object) standardGeneric("clinical"))

#' @describeIn FirehoseData Extract clinical data
setMethod("clinical", "FirehoseData", function(object) {
    getElement(object, "clinical")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("RNASeqGene", function(object) standardGeneric("RNASeqGene"))

#' @describeIn FirehoseData Extract RNASeqGene
setMethod("RNASeqGene", "FirehoseData", function(object) {
    getElement(object, "RNASeqGene")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("RNASeq2GeneNorm", function(object) standardGeneric("RNASeq2GeneNorm"))

#' @describeIn FirehoseData Extract RNASeq2GeneNorm normalized
setMethod("RNASeq2GeneNorm", "FirehoseData", function(object) {
    getElement(object, "RNASeq2GeneNorm")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("miRNASeqGene", function(object) standardGeneric("miRNASeqGene"))

#' @describeIn FirehoseData Extract miRNASeqGene microRNA
setMethod("miRNASeqGene", "FirehoseData", function(object) {
    getElement(object, "miRNASeqGene")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("CNASNP", function(object) standardGeneric("CNASNP"))

#' @describeIn FirehoseData Extract Copy Number Alterations
setMethod("CNASNP", "FirehoseData", function(object) {
    getElement(object, "CNASNP")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("CNVSNP", function(object) standardGeneric("CNVSNP"))

#' @describeIn FirehoseData Extract Copy Number Variations
setMethod("CNVSNP", "FirehoseData", function(object) {
    getElement(object, "CNVSNP")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("CNASeq", function(object) standardGeneric("CNASeq"))

#' @describeIn FirehoseData Copy Number Alterations Sequencing
setMethod("CNASeq", "FirehoseData", function(object) {
    getElement(object, "CNASeq")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("CNACGH", function(object) standardGeneric("CNACGH"))

#' @describeIn FirehoseData Copy Number Alterations Comparative Genomic
#' Hybridization
setMethod("CNACGH", "FirehoseData", function(object) {
    getElement(object, "CNACGH")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("Methylation", function(object) standardGeneric("Methylation"))

#' @describeIn FirehoseData Methylation
setMethod("Methylation", "FirehoseData", function(object) {
    getElement(object, "Methylation")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("mRNAArray", function(object) standardGeneric("mRNAArray"))

#' @describeIn FirehoseData Methylation
setMethod("mRNAArray", "FirehoseData", function(object) {
    getElement(object, "mRNAArray")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("miRNAArray", function(object) standardGeneric("miRNAArray"))

#' @describeIn FirehoseData microRNA array
setMethod("miRNAArray", "FirehoseData", function(object) {
    getElement(object, "miRNAArray")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("RPPAArray", function(object) standardGeneric("RPPAArray"))

#' @describeIn FirehoseData Reverse Phase Protein Array
setMethod("RPPAArray", "FirehoseData", function(object) {
    getElement(object, "RPPAArray")
})

#' @rdname FirehoseData-extractors
#' @export
setGeneric("Mutation", function(object) standardGeneric("Mutation"))

#' @describeIn FirehoseData Mutation
setMethod("Mutation", "FirehoseData", function(object) {
    getElement(object, "Mutation")
})

#' @rdname FirehoseData-extractors 
#' @export
setGeneric("GISTIC", function(object) standardGeneric("GISTIC"))

#' @describeIn FirehoseData GISTIC v2 thresholded by gene
setMethod("GISTIC", "FirehoseData", function(object) {
    getElement(object, "GISTIC")
})
