#' Get list of TCGA cohorts.
#'
#' \code{getFirehoseDatasets} returns a character array for cohorts.
#'
#' @return A character string
#' @examples
#' getFirehoseDatasets()
#' @export getFirehoseDatasets
getFirehoseDatasets <- function(){
  runDataset <- read.table("http://www.canevolve.org/fmineRdataset.txt",
                           header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  runDataset <- as.character(runDataset[,1])
  return(runDataset)
}

#' RTCGASample
#' 
#' Example dataset not biologically meaningful
#' @docType data
#' 
#' @keywords data
"RTCGASample"
