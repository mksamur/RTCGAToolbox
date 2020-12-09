#' Create a SummarizedExperiment from FireHose GISTIC
#'
#' @description Use the output of \code{getFirehoseData} to create a
#' \linkS4class{SummarizedExperiment}. This can be done for three types of
#' data, G-scores thresholded by gene, copy number by gene, and copy number by
#' peak regions.
#'
#' @param gistic A \link[RTCGAToolbox]{FirehoseGISTIC-class} object
#' @param dataType Either one of "ThresholdedByGene", "AllByGene", "Peaks"
#' @param ... Additional arguments passed to 'getGISTICPeaks'.
#'
#' @author L. Geistlinger, M. Ramos
#'
#' @examples
#'
#' co <- getFirehoseData("COAD", clinical = FALSE, GISTIC = TRUE,
#'     destdir = tempdir())
#' makeSummarizedExperimentFromGISTIC(co, "AllByGene")
#'
#' @return A \code{SummarizedExperiment} object
#' @export
makeSummarizedExperimentFromGISTIC <- function(gistic, dataType, ...) {
    if (identical(dataType, "Peaks")) {
        gist <- getGISTICPeaks(gistic, ...)
        if (!length(gist)) return(list())
        rel.cols <- grepl("^TCGA", colnames(gist))

        gistData <- as.matrix(gist[, rel.cols])
        # get the peak type (amplification / deletion)
        peak.type <- vapply(strsplit(gist[["Unique.Name"]], " "),
            function(x) x[[1L]], character(1L))
        rowdata <- cbind.data.frame(gist[, !rel.cols], type = peak.type,
            stringsAsFactors = FALSE)
        feats <- grepl("gene|ranges", names(rowdata), ignore.case = TRUE)
        rows <- rowdata[, feats]
        if (length(rows)) {
            if (as.logical(anyDuplicated(rows))) {
                uniq <- !duplicated(rows)
                rows <- rows[uniq]
                rowdata <- rowdata[uniq, ]
                gistData <- gistData[uniq, ]
            }
            rownames(rowdata) <- rows
        }
        rowranges <- gist[["rowRanges"]]
        rowranges <- as(rowranges, "GRanges")
        # get the peak type (amplification / deletion)
        peak.type <- vapply(strsplit(gist[["Unique.Name"]], " "),
            function(x) x[[1L]], character(1L))
        # create the SE
        gisticSE <- SummarizedExperiment(gistData,
            rowRanges = rowranges)
        rowData(gisticSE) <- rowdata
    } else {
        gist <- getData(gistic, "GISTIC", dataType)
        if (!length(gist)) return(list())
        rel.cols <- grepl("^TCGA", colnames(gist))
        gistData <- as.matrix(gist[, rel.cols])
        annoteRowDF <- gist[, !rel.cols, drop = FALSE]
        colnames(gistData) <- .stdIDs(colnames(gistData))
        gisticSE <- SummarizedExperiment(gistData, rowData = annoteRowDF)
    }
    gisticSE
}
