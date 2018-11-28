#' Download GISTIC2 peak-level data from the Firehose pipeline
#'
#' @description Access GISTIC2 level 4 copy number data through
#' \code{gdac.broadinstitute.org}
#'
#' @param dataset A TCGA cancer code
#' @param peak The peak type, select from "wide", "narrow", "full".
#' @param rm.chrX (logical default TRUE) Whether to remove observations in the
#' X chromosome
#' @param gistic2Date (character default "20160128") Data of the analysis
#' pipeline run
#' @param destdir Directory location to save the downloaded file
#' (default tempdir())
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
#' co <- getGISTICPeaks("COAD", "wide")
#' class(co)
#' head(co)[1:6]
#'
#' @export
getGISTICPeaks <- function(dataset,  peak = c("wide", "narrow", "full"),
    rm.chrX = TRUE, gistic2Date = "20160128", destdir = tempdir()) {

    BROAD.URL <- file.path("https://gdac.broadinstitute.org",
        paste("runs/analyses_", substr(gistic2Date, 1, 4),
        substr(gistic2Date, 5, 6), substr(gistic2Date, 7, 8), sep = "_"),
        "data", dataset, gistic2Date)

    GISTIC.FILE <- paste0("gdac.broadinstitute.org_", dataset, "-TP.",
        "CopyNumber_Gistic2.Level_4.", gistic2Date, "00.0.0")
    # download the tar
    url <- file.path(BROAD.URL, paste0(GISTIC.FILE, ".tar.gz"))
    tumorType <- switch(dataset, LAML = "TB", SKCM = "TM", "TP")
    url <- gsub("TP", tumorType, url, fixed = TRUE)
    basef <- file.path(GISTIC.FILE, "all_lesions.conf_99.txt")

    if (!file.exists(file.path(destdir, basef))) {
        tempFile <- tempfile(tmpdir = destdir, fileext = ".tar.gz")
        download.file(url, destfile = tempFile)
        files <- untar(tempFile, list=TRUE)
        # extract the lesions file
        basef <- files[grepl("all_lesions.conf_99.txt", files, fixed = TRUE)]
        untar(tempFile, files = basef, exdir = destdir)
        file.remove(tempFile)
    }

    # read the lesion file
    gistic <- read.delim(file.path(destdir, basef), as.is=TRUE)

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
