#' Download expression-based cancer subtypes from the Broad Institute
#'
#' @description Obtain the mRNA expression clustering results from the
#' Broad Institute for a specific cancer code (see [getFirehoseDatasets]).
#'
#' @param dataset A TCGA cancer code, e.g. "OV" for ovarian cancer
#' @param clust.alg The selected cluster algorithm, either "CNMF" or
#' "ConsensusPlus" (default "CNMF")
#'
#' @return A `data.frame` of cluster and silhouette values
#'
#' @importFrom RCurl url.exists
#' @importFrom S4Vectors isSingleString
#' @importFrom utils read.delim
#'
#' @author Ludwig Geistlinger
#'
#' @examples
#' co <- getBroadSubtypes("COAD", "CNMF")
#' head(co)
#'
#' @export
getBroadSubtypes <- function(dataset, clust.alg = c("CNMF", "ConsensusPlus"))
{
    if (!isSingleString(clust.alg))
        stop("Select a valid clustering algorithm")
    if (!isSingleString(dataset))
        stop("Enter a valid cancer code. See '?getFirehoseDatasets'")

    url <- file.path("https://gdac.broadinstitute.org/runs/analyses__latest",
        "reports/cancer", paste0(dataset, "-TP"),
        paste0("mRNA_Clustering_", clust.alg),
        paste0(dataset, "-TP.bestclus.txt"))

    rplmt <- switch(dataset, LAML = "TB", SKCM = "TM", "TP")
    url <- gsub("TP", rplmt, url)

    # check mRNA cLustering availability
    if (!RCurl::url.exists(url)) {
        warning(paste("mRNA clustering not available for",
            dataset, "- mRNAseq clustering is taken instead"))
        url <- sub("mRNA", "mRNAseq", url)
    }

    subtys <- read.delim(url, skip = 1, as.is = TRUE)
    rownames(subtys) <- subtys[, 1]
    subtys[, -1]
}
