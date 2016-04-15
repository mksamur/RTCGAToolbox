#' Extract TCGA barcode data
#'
#' This function uses the full TCGA barcode to return a data frame of the
#' data to the right of the participant identifier.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @return A dataframe with sample type, sample code, vial, portion, analyte,
#' plate, and center columns.
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @examples
#'
#' \dontrun{
#' barcodes <- colnames(getElement(a2, "RNASeqGene"))
#' right_bctable <- bcRight(barcodes)
#' }
#'
#' @export bcRight
bcRight <- function(barcodes) {
  sample_type <- as.character(
    samptab[,"Definition"][
      match(as.numeric(bcIDR(barcodes, sample = TRUE, part = FALSE)),
            samptab[,"Code"])])
  tb <- data.frame(patientids = bcIDR(barcodes),
                   sample_type,
                   sample_code = as.character(
                     bcIDR(barcodes, part=FALSE, sample = TRUE)),
                   vial = as.character(
                     substr(bcIDR(barcodes, part=FALSE, sample = TRUE,
                                  vial = TRUE), 3,3)),
                   portion = as.character(
                     substr(bcIDR(barcodes, part = FALSE, portion = TRUE),
                            1,2)),
                   analyte = as.character(
                     substr(bcIDR(barcodes, part = FALSE, portion = TRUE),
                            3,3)),
                   plate = as.character(
                     bcIDR(barcodes, part = FALSE, plate = TRUE)),
                   center = as.character(
                     bcIDR(barcodes, part = FALSE, center=TRUE)),
                   stringsAsFactors = FALSE)
  return(tb)
}