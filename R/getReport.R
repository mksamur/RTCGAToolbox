#' DEFUNCT: Draws a circle plot into working directory
#'
#' \code{getReport} draws a circle plot into your working directory to show log
#' fold changes for differentially expressed genes, copy number alterations and
#' mutations.
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param DGEResult1 Differential gene expression results object (Optional)
#' @param DGEResult2 Differential gene expression results object (Optional)
#' @param geneLocations Gene coordinates.
#'
#' @return Draws a circle plot
#' @export getReport
#' @examples
#' data(accmini)
getReport <- function(dataObject,DGEResult1=NULL,DGEResult2=NULL,geneLocations) {
  .Defunct(
    msg = "This function is no longer maintained."
  )
}
