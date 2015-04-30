#' Perform differential gene expression analysis for mRNA expression data.
#'
#' \code{getDiffExpressedGenes} returns a list that stores the results for each dataset.
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param DrawPlots A logical parameter to draw heatmaps and volcano plots.
#' @param adj.method Raw p value adjustment methods (Default "BH")
#' @param adj.pval Adjusted p value cut off for results table (Default 0.05)
#' @param raw.pval raw p value cut off for results table (Default 0.05)
#' @param logFC log fold change cut off for results table (Default 2)
#' @param hmTopUpN Max number of up regulated genes in heatmap (Default 100)
#' @param hmTopDownN Max number of down regulated genes in heatmap (Default 100)
#' @param meanFilter Mean read counts for each gene to filter not expressed genes (Default 10)
#' @return Returns a list that stores results for each dataset
#' @examples
#' data(RTCGASample)
#' dgegenes = getDiffExpressedGenes(RTCGASample)
#' dgegenes
#' showResults(dgegenes[[1]])
#' dgegenes = showResults(dgegenes[[1]])
#' head(dgegenes)
#' \dontrun{
#' }

getDiffExpressedGenes <- function(dataObject,DrawPlots=TRUE,adj.method="BH",adj.pval=0.05,raw.pval=0.05,logFC=2,
                                  hmTopUpN=100,hmTopDownN=100,meanFilter=10)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  validMatrix <- character()
  #check expression data matrices
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeq2GeneNorm)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
  if(length(dataObject@mRNAArray) > 0){validMatrix <- append(validMatrix,"mRNAArray")}
  if(length(validMatrix) == 0){stop("There is no valid expression data in the object!")}
  if(class(DrawPlots) != "logical" | is.null(DrawPlots)){stop("DrawPlots must be logical!")}
  if(is.null(adj.method) | is.na(adj.method) | (adj.method %in% c("BH","BY","holm","none"))){adj.method="BH"}
  if(is.null(adj.pval) | is.na(adj.pval) | length(adj.pval) > 1 | adj.pval > 1 | adj.pval < 0){adj.pval=0.05}
  if(is.null(raw.pval) | is.na(raw.pval) | length(raw.pval) > 1 | raw.pval > 1 | raw.pval < 0){raw.pval=0.05}
  if(is.null(logFC) | is.na(logFC) | length(logFC) > 1 | logFC < 0 ){logFC=2}
  if(is.null(hmTopUpN) | is.na(hmTopUpN) | length(hmTopUpN) > 1 | hmTopUpN < 0){hmTopUpN=100}
  if(is.null(hmTopDownN) | is.na(hmTopDownN) | length(hmTopDownN) > 1 | hmTopDownN < 0){hmTopDownN=100}
  listResults <- list()
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  for(i in validMatrix)
  {
    switch(i,
    "RNASeq"=
    {
      chkTmp <- as.numeric(dataObject@RNASeqGene[1,])
      if(all(is.wholenumber(chkTmp)) == FALSE){warning("RNASeq data does not look like raw counts! We will skip this data!")}
      else
      {
        sampleIDs <- colnames(dataObject@RNASeqGene)
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs),ncol=7))
        rownames(samplesDat) <- sampleIDs
        for(j in 1:length(sampleIDs))
        {
          tmpRow <- unlist(strsplit(sampleIDs[j],split="-"))
          samplesDat[sampleIDs[j],] <- tmpRow
        }
        sampleIDs1 <- as.character(samplesDat[,4])
        sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
        sampleIDs1 <- as.numeric(sampleIDs1)
        normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
        tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
        analysisGo <- TRUE
        if(is.null(normalSamples) | length(normalSamples) < 1)
        {
          warning("RNASeq data: There is no sample in the normal group!")
          analysisGo <- FALSE
        }
        if(is.null(tumorSamples) | length(tumorSamples) < 1)
        {
          warning("RNASeq data: There is no sample in the tumor group!")
          analysisGo <- FALSE
        }
        if(analysisGo)
        {
          meanCounts <- apply(dataObject@RNASeqGene,1,mean)
          voomMat <- dataObject@RNASeqGene[meanCounts > meanFilter,c(normalSamples,tumorSamples)]
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          v <- voom(voomMat,design,plot=DrawPlots)
          fit <- lmFit(v,design)
          cont.matrix <- makeContrasts("Tumor-Normal", levels=design)
          fit2 <- contrasts.fit(fit, cont.matrix)
          fit2 <- eBayes(fit2)
          aradeger <- topTable(fit2, adjust.method=adj.method, genelist=fit$genes, number=length(fit2))
          aradeger <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval,])
          aradeger <- aradeger[aradeger[,1] > logFC | aradeger[,1] < (-1*logFC),]
          tmpReturn <- new("DGEResult",Dataset="RNASeq",Toptable=data.frame(aradeger))
          listResults <- c(listResults,tmpReturn)
          if(DrawPlots)
          {
            volcanoplot(fit2,names=fit2$genes$ID,xlab="Log Fold Change",ylab="Log Odds",pch=16,cex=0.35)
            if(nrow(aradeger) > 2 ){
              aradeger <- aradeger[order(aradeger[,1],decreasing=TRUE),]
              if(nrow(aradeger) >= (hmTopDownN+hmTopUpN))
              {
                if(hmTopUpN > 0){topgenes <- rownames(aradeger)[1:hmTopUpN]}
                else{topgenes <- NULL}
                if(hmTopDownN > 0){bottomgenes <- rownames(aradeger)[(nrow(aradeger)- (hmTopDownN-1)):nrow(aradeger)]}
                else{bottomgenes <- NULL}
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- v[c(topgenes,bottomgenes),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- c(topgenes,bottomgenes)
                try(heatmap(v,col=bluered,scale="row",main="RNASeq",Colv=NA),silent=FALSE)
              }
              else
              {
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- v[rownames(aradeger),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- rownames(aradeger)
                try(heatmap(v,col=bluered,scale="row",main="RNASeq",Colv=NA),silent=FALSE)
              }
            }
          }
        }
      }
    },
    "RNASeq2"=
    {
      chkTmp <- as.numeric(dataObject@RNASeq2GeneNorm[1,])
      if(all(is.wholenumber(chkTmp)) == FALSE){warning("RNASeq2 data does not look like raw counts! We will skip this data!")}
      else
      {
        sampleIDs <- colnames(dataObject@RNASeq2GeneNorm)
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs),ncol=7))
        rownames(samplesDat) <- sampleIDs
        for(j in 1:length(sampleIDs))
        {
          tmpRow <- unlist(strsplit(sampleIDs[j],split="-"))
          samplesDat[sampleIDs[j],] <- tmpRow
        }
        sampleIDs1 <- as.character(samplesDat[,4])
        sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
        sampleIDs1 <- as.numeric(sampleIDs1)
        normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
        tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
        analysisGo <- TRUE
        if(is.null(normalSamples) | length(normalSamples) < 1)
        {
          warning("RNASeq2 data: There is no sample in the normal group!")
          analysisGo <- FALSE
        }
        if(is.null(tumorSamples) | length(tumorSamples) < 1)
        {
          warning("RNASeq2 data: There is no sample in the tumor group!")
          analysisGo <- FALSE
        }
        if(analysisGo)
        {
          meanCounts <- apply(dataObject@RNASeq2GeneNorm,1,mean)
          voomMat <- dataObject@RNASeq2GeneNorm[meanCounts > meanFilter,c(normalSamples,tumorSamples)]
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          v <- voom(voomMat,design,plot=DrawPlots)
          fit <- lmFit(v,design)
          cont.matrix <- makeContrasts("Tumor-Normal", levels=design)
          fit2 <- contrasts.fit(fit, cont.matrix)
          fit2 <- eBayes(fit2)
          aradeger <- topTable(fit2, adjust.method=adj.method, genelist=fit$genes, number=length(fit2))
          aradeger <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval,])
          aradeger <- aradeger[aradeger[,1] > logFC | aradeger[,1] < (-1*logFC),]
          tmpReturn <- new("DGEResult",Dataset="RNASeq2",Toptable=data.frame(aradeger))
          listResults <- c(listResults,tmpReturn)
          if(DrawPlots)
          {
            volcanoplot(fit2,names=fit2$genes$ID,xlab="Log Fold Change",ylab="Log Odds",pch=16,cex=0.35)
            if(nrow(aradeger) > 2 ){
              aradeger <- aradeger[order(aradeger[,1],decreasing=TRUE),]
              if(nrow(aradeger) >= (hmTopDownN+hmTopUpN))
              {
                if(hmTopUpN > 0){topgenes <- rownames(aradeger)[1:hmTopUpN]}
                else{topgenes <- NULL}
                if(hmTopDownN > 0){bottomgenes <- rownames(aradeger)[(nrow(aradeger)- (hmTopDownN-1)):nrow(aradeger)]}
                else{bottomgenes <- NULL}
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- v[c(topgenes,bottomgenes),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- c(topgenes,bottomgenes)
                try(heatmap(v,col=bluered,scale="row",main="RNASeq2",Colv=NA),silent=FALSE)
              }
              else
              {
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- v[rownames(aradeger),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- rownames(aradeger)
                try(heatmap(v,col=bluered,scale="row",main="RNASeq2",Colv=NA),silent=FALSE)
              }
            }
          }
        }
      }
    },
    "mRNAArray"=
    {
      for(j in 1:length(dataObject@mRNAArray))
      {
        tmpObj <- dataObject@mRNAArray[[j]]
        geneMat <- tmpObj@DataMatrix
        sampleIDs <- colnames(geneMat)
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs),ncol=7))
        rownames(samplesDat) <- sampleIDs
        for(j in 1:length(sampleIDs))
        {
          if(grepl("\\.",sampleIDs[j]))
          {
            tmpRow <- unlist(strsplit(sampleIDs[j],split="\\."))
            samplesDat[sampleIDs[j],] <- tmpRow
          }
          else
          {
            tmpRow <- unlist(strsplit(sampleIDs[j],split="-"))
            samplesDat[sampleIDs[j],] <- tmpRow
          }
        }
        sampleIDs1 <- as.character(samplesDat[,4])
        sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
        sampleIDs1 <- as.numeric(sampleIDs1)
        normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
        tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
        analysisGo <- TRUE
        if(is.null(normalSamples) | length(normalSamples) < 1)
        {
          message("mRNA array data: There is no sample in the normal group!")
          message(tmpObj@Filename)
          warning("mRNA array data: There is no sample in the normal group!")
          analysisGo <- FALSE
        }
        if(is.null(tumorSamples) | length(tumorSamples) < 1)
        {
          message("mRNA array data: There is no sample in the tumor group!")
          message(tmpObj@Filename)
          warning("mRNA array data: There is no sample in the tumor group!")
          analysisGo <- FALSE
        }
        if(analysisGo)
        {
          geneMat <- geneMat[,c(normalSamples,tumorSamples)]
          rN <- rownames(geneMat)
          cN <- colnames(geneMat)
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          fit <- lmFit(geneMat,design)
          cont.matrix <- makeContrasts("Tumor-Normal", levels=design)
          fit2 <- contrasts.fit(fit, cont.matrix)
          fit2 <- eBayes(fit2)
          aradeger <- topTable(fit2, adjust.method=adj.method, genelist=fit$genes, number=length(fit2))
          aradeger <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval,])
          aradeger <- aradeger[aradeger[,1] > logFC | aradeger[,1] < (-1*logFC),]
          tmpReturn <- new("DGEResult",Dataset=tmpObj@Filename,Toptable=data.frame(aradeger))
          listResults <- c(listResults,tmpReturn)
          if(DrawPlots)
          {
            volcanoplot(fit2,names=fit2$genes$ID,xlab="Log Fold Change",ylab="Log Odds",pch=16,cex=0.35)
            if(nrow(aradeger) > 2 ){
              aradeger <- aradeger[order(aradeger[,1],decreasing=TRUE),]
              if(nrow(aradeger) >= (hmTopDownN+hmTopUpN))
              {
                if(hmTopUpN > 0){topgenes <- rownames(aradeger)[1:hmTopUpN]}
                else{topgenes <- NULL}
                if(hmTopDownN > 0){bottomgenes <- rownames(aradeger)[(nrow(aradeger)- (hmTopDownN-1)):nrow(aradeger)]}
                else{bottomgenes <- NULL}
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- geneMat[c(topgenes,bottomgenes),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- c(topgenes,bottomgenes)
                hmHead <- paste("mRNA Array #",1,sep="")
                for(pp in 1:2)
                {
                  message("##############################")
                }
                message("Array heatmap info:")
                message(paste(hmHead,":",tmpObj@Filename))
                for(pp in 1:2)
                {
                  message("##############################")
                }
                try(heatmap(v,col=bluered,scale="row",main=hmHead,Colv=NA),silent=FALSE)
              }
              else
              {
                bluered <- colorRampPalette(c("blue","white","red"))(256)
                v <- geneMat[rownames(aradeger),]
                v <- apply(v,2,as.numeric)
                rownames(v) <- rownames(aradeger)
                try(heatmap(v,col=bluered,scale="row",main=tmpObj@Filename,Colv=NA),silent=FALSE)
              }
            }
          }
        }
      }
    })
  }
  return(listResults)
}