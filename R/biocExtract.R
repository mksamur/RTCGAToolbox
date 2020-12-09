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
#' \code{\link[RaggedExperiment]{RaggedExperiment-class}} for more details.
#'
#' @details A typical additional argument for this function passed down to
#' lower level functions is the `names.field` which indicates the row names
#' in the data. By default, it is the "Hugo_Symbol" column in the internal
#' code that converts `data.frame`s to `SummarizedExperiment` representations
#' (via the `.makeSummarizedExperimentFromDataFrame` internal function).
#'
#' @section type:
#' Choices include:
#' \itemize{
#'     \item{clinical} - Get the clinical data slot
#'     \item{RNASeqGene} - RNASeqGene - RNASeq v1
#'     \item{RNASeqGene} - RNASeq2Gene - RNASeq v2
#'     \item{RNASeq2GeneNorm} - RNASeq v2 Normalized
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
#' @param ... Additional arguments passed to lower level functions that
#' convert tabular data into Bioconductor object such as
#' \code{.makeRangedSummarizedExperimentFromDataFrame} or
#' \code{.makeRaggedExperimentFromDataFrame}
#'
#' @return Either an \linkS4class{SummarizedExperiment} object or a
#' \linkS4class{RaggedExperiment} object.
#'
#' @md
#'
#' @author Marcel Ramos \email{marcel.ramos@@roswellpark.org}
#'
#' @examples
#'
#' data(accmini)
#' biocExtract(accmini, "RNASeq2Gene")
#' biocExtract(accmini, "miRNASeqGene")
#' biocExtract(accmini, "RNASeq2GeneNorm")
#' biocExtract(accmini, "CNASNP")
#' biocExtract(accmini, "CNVSNP")
#' biocExtract(accmini, "Methylation")
#' biocExtract(accmini, "Mutation")
#' biocExtract(accmini, "RPPAArray")
#' biocExtract(accmini, "GISTIC")
#'
#' @export biocExtract
biocExtract <- function(object, type = c("clinical", "RNASeqGene",
    "RNASeq2Gene", "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP",
    "CNASeq", "CNACGH", "Methylation", "Mutation", "mRNAArray", "miRNAArray",
    "RPPAArray", "GISTIC", "GISTICA", "GISTICT", "GISTICP"), ...) {
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
    if (is.matrix(object)) {
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
        result <-
        if (type == "GISTIC") {
            names(slotreq) <- slotreq
            Filter(length,
                lapply(slotreq, function (x) {
                    makeSummarizedExperimentFromGISTIC(object, x)
                })
            )
        } else {
            makeSummarizedExperimentFromGISTIC(object, slotreq)
        }
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
        object <- .convertRangeBioc(object)
    } else {
        object <- .standardizeBC(object)
        metadat <- metadata(object)
        object <- tryCatch({
            data.matrix(object)
        }, error = function(e) {
            object
        })
        object <- SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    return(object)
}

.convertRangeBioc <- function(object, ...) {
    build <- unique(.getBuild(object))
    build <- TCGAutils::correctBuild(build, "NCBI")
    if (.hasConsistentRanges(object)) {
        object <- .makeRangedSummarizedExperimentFromDataFrame(object,
            build = build, ...)
    } else {
        split.field <- .findSampleCol(object)
        if (is.na(split.field) || !length(split.field))
            object <- .makeGRangesFromDataFrame(object,
                build = build)
        else
            object <- .makeRaggedExperimentFromDataFrame(
                object, split.field = split.field, build = build, ...)
    }
    object
}
