#' Perform correlation analysis betwwen gene expression and copy number data
#'
#' \code{getCNGECorrelation} returns a list that stores the results correlation between gene expression and copy number data.
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param adj.method Raw p value adjustment methods (Default "BH")
#' @param adj.pval Adjusted p value cut off for results table (Default 0.05)
#' @param raw.pval raw p value cut off for results table (Default 0.05)
#' @return Returns a list that stores results for each dataset
#' @examples
#' data(RTCGASample)
#' \dontrun{
#' corRes = getCNGECorrelation(a2)
#' }
getCNGECorrelation <- function(dataObject,adj.method="BH",adj.pval=0.05,raw.pval=0.05)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class! ")}
  validMatrix <- character()
  #check expression data matrices
  if(length(dataObject@mRNAArray) > 0){validMatrix <- append(validMatrix,"mRNAArray")}
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeq2GeneNorm)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
  if(is.null(adj.method) | is.na(adj.method) | (adj.method %in% c("BH","BY","holm","none"))){adj.method="BH"}
  if(is.null(adj.pval) | is.na(adj.pval) | length(adj.pval) > 1 | adj.pval > 1 | adj.pval < 0){adj.pval=0.05}
  if(is.null(raw.pval) | is.na(raw.pval) | length(raw.pval) > 1 | raw.pval > 1 | raw.pval < 0){raw.pval=0.05}
  if(length(validMatrix) == 0){stop("There is no valid expression data in the object!")}
  if(dim(dataObject@GISTIC@AllByGene)[1] == 0 | dim(dataObject@GISTIC@AllByGene)[2] == 0 ){stop("There is no GISTIC data!")}
  listResults <- list()
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  for(i in validMatrix)
  {
    if(i == "RNASeq")
    {
      chkTmp <- as.numeric(dataObject@RNASeqGene[1,])
      controlVal=FALSE
      if(all(is.wholenumber(chkTmp)) == TRUE)
      {
        #warning("Current version of correlation tool only works with normalized RNASeq data!")
        controlVal=TRUE
      }
      if(TRUE)
      {
        sampleIDs1 <- colnames(dataObject@RNASeqGene)
        sampleIDs2 <- colnames(dataObject@GISTIC@AllByGene)
        sampleIDs1 <- gsub(pattern="\\.",replacement="-",sampleIDs1)
        sampleIDs2 <- gsub(pattern="\\.",replacement="-",sampleIDs2)
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs1),ncol=7))
        rownames(samplesDat) <- sampleIDs1
        for(j in 1:length(sampleIDs1))
        {
          tmpRow <- unlist(strsplit(sampleIDs1[j],split="-"))
          samplesDat[sampleIDs1[j],] <- tmpRow
        }
        sampleIDs1 <- as.character(samplesDat[,4])
        sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
        sampleIDs1 <- as.numeric(sampleIDs1)
        samplesDat[,4] <- sampleIDs1
        sampleIDs1 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],samplesDat[,4],sep="-")
        sampleIDs2 <- sampleIDs2[4:length(sampleIDs2)]
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs2),ncol=7))
        rownames(samplesDat) <- sampleIDs2
        for(j in 1:length(sampleIDs2))
        {
          tmpRow <- unlist(strsplit(sampleIDs2[j],split="-"))
          samplesDat[sampleIDs2[j],] <- tmpRow
        }
        sampleIDs2 <- as.character(samplesDat[,4])
        sampleIDs2 <- substr(sampleIDs2,1,nchar(sampleIDs2)-1)
        sampleIDs2 <- as.numeric(sampleIDs2)
        samplesDat[,4] <- sampleIDs2
        sampleIDs2 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],samplesDat[,4],sep="-")
        tmpMat1 <- dataObject@RNASeqGene
        colnames(tmpMat1) <- sampleIDs1
        cnGenes <- dataObject@GISTIC@AllByGene[,1]
        cnGenes <- setdiff(cnGenes,cnGenes[duplicated(cnGenes)])
        tmpMat2 <- dataObject@GISTIC@AllByGene[dataObject@GISTIC@AllByGene[,1] %in% cnGenes,]
        rownames(tmpMat2) <- tmpMat2[,1]
        tmpMat2 <- tmpMat2[,4:ncol(tmpMat2)]
        colnames(tmpMat2) <- sampleIDs2
        commonSamples <- intersect(sampleIDs1,sampleIDs2)
        if(controlVal){tmpMat1=voom(tmpMat1)$E}
        if(length(commonSamples) > 5)
        {
          tmpMat1 <- tmpMat1[,commonSamples]
          tmpMat2 <- tmpMat2[,commonSamples]
          #rnaseqGenes <- rownames(tmpMat1)
          #rnaseqGenes2 <- character()
          #for(rg in rnaseqGenes)
          #{
          # rnaseqGenes2 <- append(rnaseqGenes2,as.character(strsplit(rg,"\\|")[[1]][1]))
          #}
          #rnaFrame <- data.frame(rnaseqGenes,rnaseqGenes2)
          #rnaFrame <- rnaFrame[!duplicated(rnaFrame[,2]),]
          #rnaseqGenes2 <- rnaFrame[,2]
          #names(rnaseqGenes2) <- rnaFrame[,1]
          #tmpMat1 <- tmpMat1[names(rnaseqGenes2),]
          #rownames(tmpMat1) <- rnaseqGenes2
          commonGenes <- intersect(rownames(tmpMat2),rownames(tmpMat1))
          tmpMat2 <- tmpMat2[commonGenes,]
          #message(dim(tmpMat2))
          tmpMat1 <- tmpMat1[commonGenes,]
          #message(dim(tmpMat1))
          #meanVal <- apply(tmpMat1,1,mean)
          #tmpMat1 <- tmpMat1[meanVal > summary(meanVal)[3],]
          #tmpMat2 <- tmpMat2[rownames(tmpMat1),]
          retMat <- data.frame(matrix(ncol=4,nrow=nrow(tmpMat1)))
          retMat[,1] <- as.character()
          rnaseqGenes2 <- rownames(tmpMat1)
          for(rs in 1:nrow(tmpMat1))
          {
            retMat[rs,1] <- rnaseqGenes2[rs]
            suppressWarnings(
              corTmp <- cor.test(as.numeric(tmpMat1[rs,]),as.numeric(tmpMat2[rs,]))
            )
            retMat[rs,2] <- corTmp$estimate
            retMat[rs,3] <- corTmp$p.value
          }
          pvals <- retMat[,3]
          pvalsadj <- p.adjust(pvals, method=adj.method)
          retMat[,3] <- pvalsadj
          retMat[,4] <- pvals
          colnames(retMat) <- c("GeneSymbol","Cor","adj.p.value","p.value")
          retMat <- retMat[retMat[,3] < adj.pval & retMat[,4] < raw.pval,]
          tmpReturn <- new("CorResult",Dataset="RNASeq",Correlations=retMat)
          listResults <- c(listResults,tmpReturn)
        }
      }
    }
    else if(i == "mRNAArray")
    {
      for(jj in 1:length(dataObject@mRNAArray))
      {
        sampleIDs1 <- colnames(dataObject@mRNAArray[[jj]]@DataMatrix)
        sampleIDs2 <- colnames(dataObject@GISTIC@AllByGene)
        sampleIDs1 <- gsub(pattern="\\.",replacement="-",sampleIDs1)
        sampleIDs2 <- gsub(pattern="\\.",replacement="-",sampleIDs2)
        sampleIDs1 <- sampleIDs1[1:length(sampleIDs1)]
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs1),ncol=7))
        rownames(samplesDat) <- sampleIDs1
        for(j in 1:length(sampleIDs1))
        {
          tmpRow <- unlist(strsplit(sampleIDs1[j],split="-"))
          samplesDat[sampleIDs1[j],] <- tmpRow
        }
        sampleIDs1 <- as.character(samplesDat[,4])
        sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
        sampleIDs1 <- as.numeric(sampleIDs1)
        samplesDat[,4] <- sampleIDs1
        sampleIDs1 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],samplesDat[,4],sep="-")
        sampleIDs2 <- sampleIDs2[4:length(sampleIDs2)]
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs2),ncol=7))
        rownames(samplesDat) <- sampleIDs2
        for(j in 1:length(sampleIDs2))
        {
          tmpRow <- unlist(strsplit(sampleIDs2[j],split="-"))
          samplesDat[sampleIDs2[j],] <- tmpRow
        }
        sampleIDs2 <- as.character(samplesDat[,4])
        sampleIDs2 <- substr(sampleIDs2,1,nchar(sampleIDs2)-1)
        sampleIDs2 <- as.numeric(sampleIDs2)
        samplesDat[,4] <- sampleIDs2
        sampleIDs2 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],samplesDat[,4],sep="-")
        commonSamples <- intersect(sampleIDs1,sampleIDs2)
        #commonGenes <- intersect(dataObject@mRNAArray[[jj]]@DataMatrix[,1],dataObject@GISTIC@AllByGene[,1])
        commonGenes <- intersect(rownames(dataObject@mRNAArray[[jj]]@DataMatrix),dataObject@GISTIC@AllByGene[,1])
        if(length(commonSamples) > 5)
        {
          tmpMat1 <- dataObject@mRNAArray[[jj]]@DataMatrix
          #rownames(tmpMat1) <- tmpMat1[,1]
          #tmpMat1 <- tmpMat1[,2:ncol(tmpMat1)]
          colnames(tmpMat1) <- sampleIDs1
          tmpMat1 <- tmpMat1[commonGenes,commonSamples]
          tmpMat2 <- dataObject@GISTIC@AllByGene
          rownames(tmpMat2) <- tmpMat2[,1]
          tmpMat2 <- tmpMat2[,4:ncol(tmpMat2)]
          colnames(tmpMat2) <- sampleIDs2
          tmpMat2 <- tmpMat2[commonGenes,commonSamples]
          retMat <- data.frame(matrix(ncol=4,nrow=nrow(tmpMat1)))
          retMat[,1] <- as.character()
          rnaseqGenes2 <- rownames(tmpMat2)
          for(rs in 1:nrow(tmpMat1))
          {
            retMat[rs,1] <- rnaseqGenes2[rs]
            suppressWarnings(
              corTmp <- cor.test(as.numeric(tmpMat1[rs,]),as.numeric(tmpMat2[rs,]))
            )
            retMat[rs,2] <- corTmp$estimate
            retMat[rs,3] <- corTmp$p.value
          }
          pvals <- retMat[,3]
          pvalsadj <- p.adjust(pvals, method=adj.method)
          retMat[,3] <- pvalsadj
          retMat[,4] <- pvals
          colnames(retMat) <- c("GeneSymbol","Cor","adj.p.value","p.value")
          retMat <- retMat[retMat[,3] < adj.pval & retMat[,4] < raw.pval,]
          tmpReturn <- new("CorResult",Dataset=dataObject@mRNAArray[[jj]]@Filename,Correlations=retMat)
          listResults <- c(listResults,tmpReturn)
        }
      }
    }
  }
  return(listResults)
}