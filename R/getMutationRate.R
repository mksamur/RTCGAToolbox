#' Make a table for mutation rate of each gene in the cohort
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @return Returns a data table
#' @examples
#'
#' \dontrun{
#' data(RTCGASample)
#' mutRate = getMutationRate(dataObject=a2)
#' }
getMutationRate <- function(dataObject)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData"){stop("'dataObject' must be 'FirehoseData' class!")}
  if(!is.null(dataObject@Mutations) & dim(dataObject@Mutations)[1] > 0 & dim(dataObject@Mutations)[2] > 0)
  {
    mutAll <- dataObject@Mutations
    uniqueGenes <- unique(mutAll[,1])
    uniqueSamples <- unique(mutAll[,16])
    countsMut <- matrix(0,length(uniqueGenes),length(uniqueSamples))
    rownames(countsMut) <- uniqueGenes
    colnames(countsMut) <- uniqueSamples
    for(i in uniqueSamples)
    {
      tmpMut <- as.character(mutAll[mutAll[,16]==i,1])
      countsMut[tmpMut,i] = 1
    }
    mutRatio <- rowMeans(countsMut)
    retMat <- data.frame(Genes=rownames(countsMut),MutationRatio=mutRatio)
    return(retMat)
  }
}