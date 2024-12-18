#' Create a SummarizedExperiment from FireHose GISTIC
#'
#' @description Use the output of \code{getFirehoseData} to create a
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}.
#' This can be done for three types of data, G-scores thresholded by gene, copy
#' number by gene, and copy number by peak regions.
#'
#' @param gistic A \link[RTCGAToolbox]{FirehoseGISTIC-class} object
#' @param dataType character(1) One of "ThresholdedByGene" (default),
#'   "AllByGene", or "Peaks"
#' @param rownameCol character(1) The name of the column in the data to use as
#'   rownames in the data matrix (default: 'Gene.Symbol'). The row names are
#'   only set when the column name is found in the data and all values are
#'   unique.
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
makeSummarizedExperimentFromGISTIC <-
    function(
        gistic,
        dataType = c("AllByGene", "ThresholdedByGene", "Peaks"),
        rownameCol = "Gene.Symbol", ...
    )
{
    if (!is(gistic, "FirehoseData") && !is(gistic, "FirehoseGISTIC"))
        stop("'gistic' must be a 'FirehoseData' or 'FirehoseGISTIC' object")
    dataType <- match.arg(dataType)
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
                gist <- gist[uniq, , drop = FALSE]
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
        anyduprows <- anyDuplicated(annoteRowDF[[rownameCol]])
        if (rownameCol %in% names(annoteRowDF) && !anyduprows)
            rownames(gistData) <- annoteRowDF[[rownameCol]]
        colnames(gistData) <- .stdIDs(colnames(gistData))
        gisticSE <- SummarizedExperiment(gistData, rowData = annoteRowDF)
    }
    gisticSE
}
