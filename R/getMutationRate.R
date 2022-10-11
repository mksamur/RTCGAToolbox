#' Make a table for mutation rate of each gene in the cohort
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @return Returns a data table
#' @export getMutationRate
#' @examples
#' data(accmini)
#' mutRate <- getMutationRate(dataObject=accmini)
#' mutRate <- mutRate[order(mutRate[,2],decreasing = TRUE),]
#' head(mutRate)
getMutationRate <- function(dataObject)
{
  if(is.null(dataObject) | !is(dataObject, "FirehoseData")){stop("'dataObject' must be 'FirehoseData' class!")}
  if(!is.null(dataObject@Mutation) & dim(dataObject@Mutation)[1] > 0 & dim(dataObject@Mutation)[2] > 0)
  {
    mutAll <- dataObject@Mutation
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
