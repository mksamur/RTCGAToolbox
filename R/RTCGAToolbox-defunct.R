#' @name RTCGAToolbox-defunct
#'
#' @aliases getCNGECorrelation getDiffExpressedGenes getSurvival
#'
#' @title Functions that are defunct in the package.

#' @description \code{getCNGECorrelation}: Perform correlation analysis between
#'   gene expression and copy number data. \code{getDiffExpressedGenes}: Perform
#'   differential gene expression analysis for mRNA expression data.
#'   \code{getSurvival}: Perform survival analysis based on gene expression
#'   data.
#'
#' @details
#' \code{getCNGECorrelation} returns a list that stores the results correlation
#' between gene expression and copy number data.
#' \code{getDiffExpressedGenes} returns a list that stores the results for each
#' dataset.
#' \code{getSurvival} draws a KM plot and show survival analysis results between
#' groups that are defined by gene expression data
#'
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param adj.method Raw p value adjustment methods (Default "BH")
#' @param adj.pval Adjusted p value cut off for results table (Default 0.05)
#' @param raw.pval raw p value cut off for results table (Default 0.05)
#' @return Returns a list that stores results for each dataset
#' @export
getCNGECorrelation <- function(dataObject,adj.method="BH",adj.pval=0.05,raw.pval=0.05)
{
  .Defunct(
    msg = "This function is no longer maintained and is defunct."
  )
}

#' @rdname RTCGAToolbox-defunct
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param DrawPlots A logical parameter to draw heatmaps and volcano plots.
#' @param adj.method Raw p value adjustment methods (Default "BH")
#' @param adj.pval Adjusted p value cut off for results table (Default 0.05)
#' @param raw.pval raw p value cut off for results table (Default 0.05)
#' @param logFC log fold change cut off for results table (Default 2)
#' @param hmTopUpN Max number of up regulated genes in heatmap (Default 100)
#' @param hmTopDownN Max number of down regulated genes in heatmap (Default 100)
#' @param meanFilter Mean read counts for each gene to filter not expressed
#'   genes (Default 10)
#'
#' @return getDiffExpressedGenes: Returns a list that stores results for each
#'   dataset
#'
#' @export
getDiffExpressedGenes <- function(
    dataObject,DrawPlots=TRUE,adj.method="BH",adj.pval=0.05,raw.pval=0.05,
    logFC=2, hmTopUpN=100,hmTopDownN=100,meanFilter=10)
{
  .Defunct(
    msg = "This function is no longer maintained and is defunct."
  )
}

#' @rdname RTCGAToolbox-defunct
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param numberofGroups Can be set as 2 or 3. (Default 2) Order and divide samples into n groups by using gene expression data.
#' @param geneSymbols Gene symbol that is going to be tested
#' @param sampleTimeCensor a data frame that stores clinical data. First column should store sample IDs, second column should have time and third column should have event information. For more information please see vignette.
#' @return getSurvival: Draws a KM plot
#' @export
getSurvival <- function(dataObject,numberofGroups=2,geneSymbols,sampleTimeCensor)
{
  .Defunct(
    msg = "This function is no longer maintained and is defunct."
  )
}
