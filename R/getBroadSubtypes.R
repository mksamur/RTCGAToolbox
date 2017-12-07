#' Download Broad Institute Subtype Clustering Results
#'
#' @description Obtain the subtype clustering results from the Broad for
#' a specific cancer code (see \link{getFirehoseDatasets}).
#'
#' @param dataset A TCGA cancer code
#' @param clust.alg The selected cluster algorithm, either "CNMF" or
#' "Consensu_Plus"
#'
#' @return A \code{data.frame} of cluster and silhouette values
#'
#' @importFrom RCurl url.exists
#'
#' @author Ludwig Geistlinger
#'
#' @examples
#' co <- getBroadSubtypes("COAD", "CNMF")
#' head(co)
#'
#' @export
getBroadSubtypes <- function(dataset, clust.alg = c("CNMF", "Consensus_Plus"))
{
    stopifnot(S4Vectors::isSingleString(clust.alg), S4Vectors::isSingleString(dataset))

    url <- file.path("http://gdac.broadinstitute.org/runs/analyses__latest",
        "reports/cancer", paste0(dataset, "-TP"),
        paste0("mRNA_Clustering_", clust.alg),
        paste0(dataset, "-TP.bestclus.txt"))

    if (dataset == "LAML") url <- gsub("TP", "TB", url)
    else if (dataset == "SKCM") url <- gsub("TP", "TM", url)

    # check mRNA cLustering availability
    if (!RCurl::url.exists(url)) {
        warning(paste("mRNA clustering not available for",
            dataset, "- mRNAseq clustering is taken instead"))
        url <- sub("mRNA", "mRNAseq", url)
    }

    subtys <- read.delim(url, skip=1, as.is=TRUE)
    subtys[,1] <- TCGAutils::TCGAbarcode(subtys[, 1], sample=TRUE)
    subtys[,1] <- sub("A$", "", subtys[,1])
    rownames(subtys) <- subtys[,1]
    subtys[,-1]
}
