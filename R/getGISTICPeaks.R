#' Download peak-level data from the Firehose pipeline
#'
#' @description Access copy number GISTIC2 level 4 data through
#' \code{gdac.broadinstitute.org}
#'
#' @param dataset A TCGA cancer code
#' @param peak The peak type, select from "wide", "narrow", "full".
#' @param rm.chrX (logical default TRUE) Whether to remove observations in the
#' X chromosome
#' @param gistic2Date (character default "20160128") Data of the analysis
#' pipeline run
#' @importFrom GenomicRanges seqnames order
#' @importFrom GenomeInfoDb orderSeqlevels seqlevels seqlevels<-
#'
#' @author Ludwig Geistlinger
#'
#' @examples
#'
#' co <- getGISTICPeaks("COAD", "wide")
#' class(co)
#' head(co)[1:6]
#'
#' @export
getGISTICPeaks <- function(dataset,  peak = c("wide", "narrow", "full"),
    rm.chrX = TRUE, gistic2Date = "20160128") {

    BROAD.URL <- file.path("https://gdac.broadinstitute.org",
        paste("runs/analyses_", substr(gistic2Date, 1, 4),
        substr(gistic2Date, 5, 6), substr(gistic2Date, 7, 8), sep = "_"),
        "data", dataset, gistic2Date)

    GISTIC.FILE <- paste0("gdac.broadinstitute.org_", dataset, "-TP.",
        "CopyNumber_Gistic2.Level_4.", gistic2Date, "00.0.0.tar.gz")

    # download the tar
    url <- file.path(BROAD.URL, GISTIC.FILE)
    tumorType <- switch(dataset, LAML = "TB", SKCM = "TM", "TP")
    url <- gsub("TP", tumorType, url, fixed = TRUE)

    tempFile <- tempfile(fileext = ".tar.gz")
    download.file(url, destfile = tempFile)
    files <- untar(tempFile, list=TRUE)

    # extract the lesions file
    basef <- files[grepl("all_lesions.conf_99.txt", files, fixed = TRUE)]
    untar(tempFile, files = basef, exdir = dirname(tempFile))

    # read the lesion file
    gistic <- read.delim(file.path(dirname(tempFile), basef), as.is=TRUE)
    unlink(dirname(tempFile), recursive=TRUE, force=TRUE)

    validCols <- vapply(gistic, function(x) !all(is.na(x)), logical(1L))
    rel.rows <- grepl("Peak +[0-9]+$", gistic[["Unique.Name"]])
    gistic <- gistic[rel.rows, validCols]

    isBCodes <- grepl("^TCGA", names(gistic))
    names(gistic)[isBCodes] <- .stdIDs(names(gistic)[isBCodes])

    # (a) get the ranges from chosen peaks
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
