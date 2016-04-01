#' Convert raw TCGA Mutation data to GRangesList
#' 
#' This function takes the data.frame of raw data from the output of a TCGA
#' data pipeline and converts it to a \linkS4class{GRangesList} class.
#' 
#' @param inputData A \code{data.frame} class of TCGA mutation data.
#' @return A \linkS4class{GRangesList} class object
#' 
#' @author Marcel Ramos \email{mramos09@gmail.com}
#' 
#' @export
makeGRangesList <- function(inputData) {
    if (!is(inputData, "data.frame")) {
        stop("inputData is not a data.frame")
    }
    names(inputData) <- tolower(names(inputData))
    longNames <- c("chromosome", "start_position", "end_position")
    shortNames <- c("chrom", "start", "end")
    twoMeta <- ifelse(all(c("num_probes", "segment_mean") %in%
                              names(inputData)), TRUE, FALSE)
    hugo <- ifelse("hugo_symbol" %in% names(inputData), TRUE, FALSE)
    if ("ncbi_build" %in% names(inputData)) {
        ncbi_build <- Reduce(intersect, inputData$ncbi_build)
        if (length(ncbi_build) == 1L) {
            ncbi <- ncbi_build
        } else {
            message("NCBI build was not consistent")
        }
    }
    names(inputData) <- plyr::mapvalues(names(inputData),
                                        longNames, shortNames,
                                        warn_missing = FALSE)
    if (!all(grepl("chr", inputData$chrom[1:5], ignore.case = TRUE))) {
        inputData$chrom <- paste0("chr", inputData$chrom)
    }
    sampleIndicator <- ifelse(is.null(inputData$sample),
                              "tumor_sample_barcode", "sample")
    ## Convert data to list for GRangesList
    inputData <- split(inputData, bcIDR(as.character(
        inputData[, sampleIndicator]),
        sample = TRUE, collapse = TRUE))
    metadats <- lapply(inputData, FUN = function(mydata) {
        mydata <- mydata[, !names(mydata)
                         %in%
                             c("seqnames", "ranges", "strand", "seqlevels",
                               "seqlengths", "isCircular", "start", "end",
                               "width", "element", "dataset", "chrom",
                               "center", "gene", "type", "chr",
                               "ref_allele", "tum_allele1", "tum_allele2")]
        return(mydata)
    })
    mygrl <- GRangesList(lapply(inputData, FUN = function(gr){
        NewGR <- GRanges(gr[, shortNames])
        if (twoMeta) {
            mcols(NewGR)$num_probes <- gr[, "num_probes"]
            mcols(NewGR)$segment_mean <- gr[, "segment_mean"]
        }
        if (hugo) { names(NewGR) <- gr[, "hugo_symbol"] }
        return(NewGR)
    }
    ))
    if (exists("ncbi")) {
        genome(mygrl) <- ncbi
    }
    metadata(mygrl) <- metadats
    return(mygrl)
}
