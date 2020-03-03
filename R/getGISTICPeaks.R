#' Download GISTIC2 peak-level data from the Firehose pipeline
#'
#' @description Access GISTIC2 level 4 copy number data through
#' \code{gdac.broadinstitute.org}
#'
#' @param object A FirehoseData GISTIC type object
#'
#' @param peak The peak type, select from "wide", "narrow", "full".
#'
#' @param rm.chrX (logical default TRUE) Whether to remove observations in the
#' X chromosome
#'
#' @return A \code{data.frame} of peak values
#'
#' @importFrom GenomicRanges seqnames order
#' @importFrom GenomeInfoDb orderSeqlevels seqlevels seqlevels<-
#' @importFrom methods as
#'
#' @author Ludwig Geistlinger
#'
#' @examples
#'
#' co <- getFirehoseData("COAD", clinical = FALSE, GISTIC = TRUE)
#' peaks <- getGISTICPeaks(co, "wide")
#' class(peaks)
#' head(peaks)[1:6]
#'
#' @export
getGISTICPeaks <-
    function(
        object,  peak = c("wide", "narrow", "full"), rm.chrX = TRUE
    )
{
    stopifnot(is(object, "FirehoseData"), is(object@GISTIC, "FirehoseGISTIC"))

    gistic <- getData(object, "GISTIC", "Peaks")

    validCols <- vapply(gistic, function(x) !all(is.na(x)), logical(1L))
    rel.rows <- grepl("Peak +[0-9]+$", gistic[["Unique.Name"]])
    gistic <- gistic[rel.rows, validCols]

    isBCodes <- grepl("^TCGA", names(gistic))
    names(gistic)[isBCodes] <- .stdIDs(names(gistic)[isBCodes])

    # (a) get the ranges from chosen peaks
    peak <- match.arg(peak, c("wide", "narrow", "full"))
    peak.col <- switch (peak,
        wide = "Wide.Peak.Limits",
        narrow = "Peak.Limits",
        full = "Region.Limits"
    )
    ranges <- gistic[, peak.col]
    ranges <- sub("\\(probes [0-9]+:[0-9]+\\) *$", "", ranges)
    grs <- as(ranges, "GRanges")

    if (rm.chrX) {
    chrx_ind <- seqnames(grs) != "chrX"
    grs <- grs[chrx_ind]
    gistic <- gistic[as.logical(chrx_ind),]
    ranges <- ranges[as.logical(chrx_ind)]
    }

    gistic <- cbind.data.frame(rowRanges = ranges, gistic,
        stringsAsFactors = FALSE)

    seq_ind <- orderSeqlevels(seqlevels(grs))
    seqlevels(grs) <- seqlevels(grs)[seq_ind]
    range_ind <- order(grs)

    gistic[range_ind, ]
}
