#' @include utils.R
NULL

#' Extract and convert data from a \code{FirehoseData} object to a
#' \code{Bioconductor} object
#'
#' This function processes data from a \code{\linkS4class{FirehoseData}} object.
#' Raw data is converted to a conventional Bioconductor object. The function
#' returns either a \linkS4class{SummarizedExperiment} or a
#' \linkS4class{RaggedExperiment} class object. In cases where there are
#' multiple platforms in a data type, an attempt to consolidate datasets will
#' be made based on matching dimension names. For ranged data, this
#' functionality is provided with more control as part of the
#' \code{RaggedExperiment} features. See
#' \code{\link[RaggedExperiment]{RaggedExperiment}} for more details.
#'
#' @section type:
#' Choices include:
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
#'     \item{GISTICA} - GISTIC v2 ('AllByGene' only)
#'     \item{GISTICT} - GISTIC v2 ('ThresholdedByGene' only)
#'     \item{GISTICP} - GISTIC v2 ('Peaks' only)
#'     \item{GISTIC} - GISTIC v2 scores, probabilities, and peaks
#' }
#'
#' @param object A \code{FirehoseData} object from which to extract data.
#' @param type The type of data to extract from the "FirehoseData" object,
#' see type section.
#' @return Either an \linkS4class{SummarizedExperiment} object or a
#' \linkS4class{RaggedExperiment} object.
#'
#' @author Marcel Ramos \email{marcel.ramos@@roswellpark.org}
#'
#' @examples \dontrun{
#'     coadmut <- getFirehoseData("COAD", Mutation = TRUE)
#'     biocExtract(coadmut, "Mutation")
#' }
#' @export biocExtract
biocExtract <- function(object, type = c("clinical", "RNASeqGene",
    "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP", "CNASeq",
    "CNACGH", "Methylation", "Mutation", "mRNAArray", "miRNAArray",
    "RPPAArray", "GISTIC", "GISTICA", "GISTICT", "GISTICP")) {
    if (!identical(length(type), 1L))
        stop("Please specify a single data type")
    message("working on: ", type)
    if (!is(object, "DataFrame") && !is.data.frame(object))
        object <- .removeShell(object, type)
    if (!length(object)) { return(object) }
    if (is.list(object) && !is.data.frame(object)) {
        object <- .unNestList(object)
    }
    if (type == "clinical") { return(object) }
    if (is(object, "matrix")) {
        return(SummarizedExperiment(assays = SimpleList(object)))
    }
    if (is(object, "list") && !is(object, "DataFrame") &&
        type != "Methylation") {
        return(.extractList(object, type = type))
    }
    if (is(object, "SummarizedExperiment")) { return(object) }

    gisticType <- grepl("^GISTIC", type, ignore.case = TRUE)
    if (gisticType) {
        slotreq <- switch(type,
            GISTICA = "AllByGene", GISTICT = "ThresholdedByGene",
            GISTICP = "Peaks",
            GISTIC = c("AllByGene", "ThresholdedByGene", "Peaks")
        )
        result <- if (type == "GISTIC") {
            names(slotreq) <- slotreq
            Filter(length,
                lapply(slotreq, function (x) { .getGISTIC(object, x) }))
        } else { .getGISTIC(object, slotreq) }
        if (length(result) == 1L)
            result <- result[[1L]]
        return(result)
    }

    if (type == "Methylation") {
        if (is.list(object)) {
            result <- lapply(object, .getMethyl)
        } else {
            result <- .getMethyl(object)
        }
        return(result)
    }

    hasRanged <- .hasRangeNames(object)
    if (hasRanged) {
        if (.hasInfo(object, "NCBI_Build")) {
            GBuild <- .getBuild(object)
        }
        if (.hasConsistentRanges(object)) {
            object <- .makeRangedSummarizedExperimentFromDataFrame(object,
                build = if (exists("GBuild")) { GBuild } else { NA })
        } else {
            object <- .makeRaggedExperimentFromDataFrame(object,
                build = if (exists("GBuild")) { GBuild } else { NA })
        }
    } else {
        object <- .standardizeBC(object)
        metadat <- metadata(object)
        object <- SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    return(object)
}
