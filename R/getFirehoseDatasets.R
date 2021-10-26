#' Get a list of TCGA disease cohorts
#'
#' \code{getFirehoseDatasets} returns a character vector of TCGA disease codes.
#' A reference table can be seen at \url{https://gdac.broadinstitute.org/}.
#'
#' @seealso \url{https://gdac.broadinstitute.org/}
#'
#' @return A character string
#'
#' @examples
#'
#' getFirehoseDatasets()
#'
#' @export getFirehoseDatasets
getFirehoseDatasets <- function() {
    c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COADREAD", "COAD",
    "DLBC", "ESCA", "FPPP", "GBMLGG", "GBM", "HNSC", "KICH", "KIPAN",
    "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO",
    "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "STES", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
}

