#' Draws a circle plot into working directory
#'
#' \code{getReport} draws a circle plot into your workin director to show log fold changes for differentially expressed genes, copy number alterations and mutations.
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param DGEResult1 Differential gene expression results object (Optional)
#' @param DGEResult2 Differential gene expression results object (Optional)
#' @param geneLocations Gene coordinates.
#' @return Draws a circle plot
#' @examples
#' data(RTCGASample)
#' data(hg19.ucsc.gene.locations)
#' t1=getDiffExpressedGenes(a2)
#' \dontrun{
#' getReport(dataObject=a2,DGEResult1=t1[[1]],geneLocations=hg19.ucsc.gene.locations)
#' }
getReport <- function(dataObject,DGEResult1=NULL,DGEResult2=NULL,geneLocations)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  if(!is.null(DGEResult1) & class(DGEResult1) != "DGEResult"){stop("DGEResult1 must be DGEResult class!")}
  if(!is.null(DGEResult2) & class(DGEResult2) != "DGEResult"){stop("DGEResult2 must be DGEResult class!")}
  pdf(file=paste(dataObject@Dataset,"-reportImage.pdf",sep=""),height=30,width=30)
  plotpos = 1;
  #require("RCircos")
  data(UCSC.HG19.Human.CytoBandIdeogram,package = "RCircos")
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL, 3, 3);
  params <- RCircos.Get.Plot.Parameters();
  params$radius.len <- 3.0;
  params$track.background <- "white"
  params$track.height <- 0.4
  params$point.size <- 2
  params$text.size <- 3
  params$track.out.start <- 0.05
  RCircos.Reset.Plot.Parameters(params)
  RCircos.Set.Plot.Area();
  RCircos.Chromosome.Ideogram.Plot();
  if(!is.null(DGEResult1))
  {
    #if(DGEResult1@Dataset=="RNASeq" | DGEResult1@Dataset=="RNASeq2")
    #{
    rnaseqGenes <- rownames(DGEResult1@Toptable)
    rnaseqGenes2 <- character()
    for(rg in rnaseqGenes)
    {
      rnaseqGenes2 <- append(rnaseqGenes2,as.character(strsplit(rg,"\\|")[[1]][1]))
    }
    rnaFrame <- data.frame(rnaseqGenes,rnaseqGenes2)
    rnaFrame <- rnaFrame[!duplicated(rnaFrame[,2]),]
    rownames(rnaFrame) <- rnaFrame[,2]
    intGenes <- intersect(rnaFrame[,2],rownames(geneLocations))
    rnaFrame <- rnaFrame[intGenes,]
    if(length(intGenes > 0))
    {
      logFC <- as.numeric(DGEResult1@Toptable[as.character(rnaFrame[,1]),1])
      histData <- cbind(geneLocations[intGenes,c(2,4,5)],logFC)
      RCircos.Scatter.Plot(histData, data.col=4, track.num=plotpos, side="in",by.fold=0.0001);
      message(paste("Track No:",plotpos," (in) differential gene expression data 1"))
      plotpos = plotpos + 1
    }
    #}
  }
  if(!is.null(DGEResult2))
  {
    #if(DGEResult1@Dataset=="RNASeq" | DGEResult1@Dataset=="RNASeq2")
    #{
    rnaseqGenes <- rownames(DGEResult1@Toptable)
    rnaseqGenes2 <- character()
    for(rg in rnaseqGenes)
    {
      rnaseqGenes2 <- append(rnaseqGenes2,as.character(strsplit(rg,"\\|")[[1]][1]))
    }
    rnaFrame <- data.frame(rnaseqGenes,rnaseqGenes2)
    rnaFrame <- rnaFrame[!duplicated(rnaFrame[,2]),]
    rownames(rnaFrame) <- rnaFrame[,2]
    intGenes <- intersect(rnaFrame[,2],rownames(geneLocations))
    rnaFrame <- rnaFrame[intGenes,]
    if(length(intGenes > 0))
    {
      logFC <- as.numeric(DGEResult1@Toptable[as.character(rnaFrame[,1]),1])
      histData <- cbind(geneLocations[intGenes,c(2,4,5)],logFC)
      RCircos.Scatter.Plot(histData, data.col=4, track.num=plotpos, side="in",by.fold=0.0001);
      message(paste("Track No:",plotpos," (in) differential gene expression data 2"))
      plotpos = plotpos + 1
    }
    #}
  }
  if(!is.null(dataObject@GISTIC) & class(dataObject@GISTIC)=="FirehoseGISTIC" & length(dataObject@GISTIC@Dataset) > 0)
  {
    cnMat <- dataObject@GISTIC@ThresholdedByGene
    rownames(cnMat) <- cnMat[,1]
    cnMat <- cnMat[,4:ncol(cnMat)]
    intGenes <- intersect(rownames(cnMat),rownames(geneLocations))
    if(length(intGenes > 0))
    {
      cnMat <- cnMat[intGenes,]
      cnMat <- apply(cnMat,2,as.numeric)
      rownames(cnMat) <- intGenes
      cnMat2 <- rowMeans(cnMat)
      cnMat <- cnMat[cnMat2 > 0.3 | cnMat2 < -0.3, ]
      #message(length(cnMat2[cnMat2 > 0.3 | cnMat2 < -0.3]))
      histData <- cbind(geneLocations[rownames(cnMat),c(2,4,5,1)],cnMat2[cnMat2 > 0.3 | cnMat2 < -0.3])
      RCircos.Heatmap.Plot(histData, data.col=5, track.num=plotpos, side="in");
      message(paste("Track No:",plotpos," (in) copy number data"))
      plotpos = plotpos + 1
    }
  }
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
    mutGenes <- rownames(countsMut)[mutRatio > 0.05]
    intGenes <- intersect(mutGenes,rownames(geneLocations))
    #message(length(intGenes))
    if(length(intGenes > 0))
    {
      histData <- geneLocations[intGenes,c(2,4,5,1)]
      RCircos.Gene.Connector.Plot(histData, track.num=1, side="out");
      RCircos.Gene.Name.Plot(histData, track.num=2,name.col=4, side="out");
      message(paste("Outside track mutations!"))
    }
  }
  dev.off()
}
