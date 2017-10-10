#' @include utils.R
NULL

#' Extract data from \code{FirehoseData} object into \code{ExpressionSet} or
#' \code{GRangesList} object
#'
#' This function processes data from a \code{\link{FirehoseData}} object.
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
#' Choices include: "RNAseqGene", "clinical", "miRNASeqGene",
#' "RNASeq2GeneNorm", "CNASNP", "CNVSNP", "CNASeq", "CNACGH", "Methylation",
#' "Mutation", "mRNAArray", "miRNAArray", "RPPAArray", "GISTIC", "GISTICA",
#' "GISTICT".
#' The "GISTICA" type of dataset represents GISTIC data by all
#' genes. "GISTICT" represents data thresholded by genes. To get both types
#' in a list, use "GISTIC".
#'
#' @param object A \code{FirehoseData} object from which to extract data.
#' @param type The type of data to extract from the "FirehoseData" object,
#' see type section.
#' @return Either an \linkS4class{SummarizedExperiment} object or a
#' \linkS4class{RaggedExperiment} object.
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @examples \dontrun{
#' library(RTCGAToolbox)
#' dataFolder <- normalizePath("~/Documents/data")
#' coadmut <- getFirehoseData("COAD", runDate = "20151101", Mutation = TRUE,
#'                          destdir = dataFolder)
#' cm <- biocExtract(coadmut, "Mutation")
#' }
#' @export biocExtract
biocExtract <- function(object, type = c("Clinical", "RNASeqGene",
    "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP", "CNASeq",
    "CNACGH", "Methylation", "Mutation", "mRNAArray", "miRNAArray",
    "RPPAArray", "GISTIC", "GISTICA", "GISTICT"), ...) {
    if (length(type) != 1L)
        stop("Please specify a single data type")
    message("working on: ", type)
    if (!is(object, "DataFrame") && !is.data.frame(object))
        object <- .removeShell(object, type)
    if (!length(object)) { return(object) }
    if (is.list(object) && !is.data.frame(object)) {
        object <- .unNestList(object)
    }
    if (type == "Clinical") { return(object) }
    if (is(object, "matrix")) {
        return(SummarizedExperiment(assays = SimpleList(object)))
    }
    if (is(object, "list") && !is(object, "DataFrame") &&
        type != "Methylation") {
        return(.extractList(object, type = type, ...))
    }
    if (is(object, "SummarizedExperiment")) { return(object) }

    gisticType <- grepl("^GISTIC", type, ignore.case = TRUE)
    if (gisticType) {
        slotreq <- switch(type, GISTICA = "AllByGene",
                          GISTICT = "ThresholdedByGene",
                          GISTIC = c("AllByGene", "ThresholdedByGene"))
        if (type == "GISTIC") {
            names(slotreq) <- slotreq
            result <- lapply(slotreq, function (x) { .getGISTIC(object, x) })
        } else
            result <- .getGISTIC(object, slotreq)
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
        if (.hasBuildInfo(object)) {
            GBuild <- .getBuild(object)
        }
        if (.hasConsistentRanges(object)) {
            object <- .makeRangedSummarizedExperimentFromDataFrame(object,
                build = if (exists("GBuild")) { GBuild } else { NULL })
        } else {
            object <- .makeRaggedExperimentFromDataFrame(object,
                build = if (exists("GBuild")) { GBuild } else { NULL })
        }
    } else {
        object <- .standardizeBC(object)
        metadat <- metadata(object)
        object <- SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    return(object)
}
