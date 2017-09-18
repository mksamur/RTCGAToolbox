#' Perform survival analysis based on gene expression data
#'
#' \code{getSurvival} draws a KM plot and show survival analysis results between groups that are defined by gene expression data
#'
#' @param dataObject This must be \code{FirehoseData} object.
#' @param numberofGroups Can be set as 2 or 3. (Default 2) Order and divide samples into n groups by using gene expression data.
#' @param geneSymbols Gene symbol that is going to be tested
#' @param sampleTimeCensor a data frame that stores clinical data. First column should store sample IDs, second column should have time and third column should have event information. For more information please see vignette.
#' @return Draws a KM plot
#' @examples
#' ## get data with  getFirehoseData() function and call survival analysis
#' ## Always check clinical data file for structural changes
#' data(RTCGASample)
#' clinicData <- getData(RTCGASample,"clinical")
#' clinicData = clinicData[,3:5]
#' clinicData[is.na(clinicData[,3]),3] = clinicData[is.na(clinicData[,3]),2]
#' survData <- data.frame(Samples=rownames(clinicData),Time=as.numeric(clinicData[,3]),
#' Censor=as.numeric(clinicData[,1]))
#' getSurvival(dataObject=RTCGASample,geneSymbols=c("FCGBP"),sampleTimeCensor=survData)
getSurvival <- function(dataObject,numberofGroups=2,geneSymbols,sampleTimeCensor)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  validMatrix <- character()
  #check expression data matrices
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeq2GeneNorm)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
  if(length(dataObject@mRNAArray) > 0){validMatrix <- append(validMatrix,"mRNAArray")}
  if(length(validMatrix) == 0){stop("There is no valid expression data in the object!")}
  if(class(numberofGroups)!="numeric"){stop("numberofGroups must be numeric!")}
  if(as.integer(numberofGroups) < 2 | as.integer(numberofGroups) > 3){stop("numberofGroups must be 2 or 3!")}
  stcs <- as.character(sampleTimeCensor[,1])
  stcs <- gsub(pattern="\\.",replacement="-",stcs)
  stcs <- toupper(stcs)
  samplesDat <- data.frame(matrix(nrow=length(stcs),ncol=3))
  rownames(samplesDat) <- stcs
  for(j in 1:length(stcs))
  {
    tmpRow <- unlist(strsplit(stcs[j],split="-"))
    samplesDat[stcs[j],] <- tmpRow
  }
  stcs <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
  rownames(sampleTimeCensor) <- stcs
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  for(i in validMatrix)
  {
    switch(i,
    "RNASeq"=
    {
      chkTmp <- as.numeric(dataObject@RNASeqGene[1,])
      controlVal = FALSE
      if(all(is.wholenumber(chkTmp)) == TRUE)
      {
        #warning("Current version of survival tool only works with normalized RNASeq data!")
        controlVal=TRUE
      }
      tmpMat1 <- dataObject@RNASeqGene
      sampleIDs1 <- colnames(tmpMat1)
      sampleIDs1 <- gsub(pattern="\\.",replacement="-",sampleIDs1)
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
      sampleIDs11 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
      colnames(tmpMat1) <- sampleIDs1
      tmpMat1 <- tmpMat1[,!duplicated(sampleIDs11)]
      colnames(tmpMat1) <- sampleIDs11[!duplicated(sampleIDs11)]
      if(controlVal){tmpMat1=voom(tmpMat1)$E}
      for(myG in geneSymbols)
      {
        if(!is.na(any(match(rownames(tmpMat1),myG))) & any(match(rownames(tmpMat1),myG)))
        {
          tmpGeneExp <- as.numeric(tmpMat1[myG,])
          names(tmpGeneExp) <- colnames(tmpMat1)
          commonSamples <- intersect(rownames(sampleTimeCensor),names(tmpGeneExp))
          if(length(commonSamples) > 9)
          {
            if(numberofGroups==2)
            {
              g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[3]]
              g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[3]]
              time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),2])
              censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),3])
              surv.group <- rep (1:2, c(length(g1s), length(g2s)))
              surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
              surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
              pvalue <- 1- pchisq (surv.diff$chisq[1], df=1)
              plot(surv.fit, xlab="Time", ylab="Survival", main=paste("RNASeq -",myG), col=c(2,4))
              pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
              legTxt <- c("< Median", ">= Median", pValueTxt)
              legend("topright", legend=legTxt, col=c(2,4,0), lty=c(1,1),cex=0.7)
            }
            else if(numberofGroups==3)
            {
              g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[2]]
              g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[2] & tmpGeneExp <= summary(tmpGeneExp)[5]]
              g3s <- names(tmpGeneExp)[tmpGeneExp > summary(tmpGeneExp)[5]]
              time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),2])
              censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),3])
              surv.group <- rep (1:3, c(length(g1s), length(g2s), length(g3s)))
              surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
              surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
              pvalue <- 1- pchisq (surv.diff$chisq[1], df=2)
              plot(surv.fit, xlab="Time", ylab="Survival", main=paste("RNASeq -",myG), col=c(2,3,4))
              pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
              legTxt <- c("< 1st Q.", "Between 1st&3rd Q.", "> 3rd. Q.", pValueTxt)
              legend("topright", legend=legTxt, col=c(2,3,4,0), lty=c(1,1),cex=0.7)
            }
          }
        }
      }
      #}
    },
    "RNASeq2"=
    {
      chkTmp <- as.numeric(dataObject@RNASeq2GeneNorm[1,])
      controlVal = FALSE
      if(all(is.wholenumber(chkTmp)) == TRUE)
      {
        #warning("Current version of survival tool only works with normalized RNASeq data!")
        controlVal=TRUE
      }
      tmpMat1 <- dataObject@RNASeq2GeneNorm
      sampleIDs1 <- colnames(tmpMat1)
      sampleIDs1 <- gsub(pattern="\\.",replacement="-",sampleIDs1)
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
      sampleIDs11 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
      colnames(tmpMat1) <- sampleIDs1
      tmpMat1 <- tmpMat1[,!duplicated(sampleIDs11)]
      colnames(tmpMat1) <- sampleIDs11[!duplicated(sampleIDs11)]
      if(controlVal){tmpMat1=voom(tmpMat1)$E}
      for(myG in geneSymbols)
      {
        if(!is.na(any(match(rownames(tmpMat1),myG))) & any(match(rownames(tmpMat1),myG)))
        {
          tmpGeneExp <- as.numeric(tmpMat1[myG,])
          names(tmpGeneExp) <- colnames(tmpMat1)
          commonSamples <- intersect(rownames(sampleTimeCensor),names(tmpGeneExp))
          if(length(commonSamples) > 9)
          {
            if(numberofGroups==2)
            {
              g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[3]]
              g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[3]]
              time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),2])
              censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),3])
              surv.group <- rep (1:2, c(length(g1s), length(g2s)))
              surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
              surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
              pvalue <- 1- pchisq (surv.diff$chisq[1], df=1)
              plot(surv.fit, xlab="Time", ylab="Survival", main=paste("RNASeq2 -",myG), col=c(2,4))
              pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
              legTxt <- c("< Median", ">= Median", pValueTxt)
              legend("topright", legend=legTxt, col=c(2,4,0), lty=c(1,1),cex=0.7)
            }
            else if(numberofGroups==3)
            {
              g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[2]]
              g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[2] & tmpGeneExp <= summary(tmpGeneExp)[5]]
              g3s <- names(tmpGeneExp)[tmpGeneExp > summary(tmpGeneExp)[5]]
              time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),2])
              censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),3])
              surv.group <- rep (1:3, c(length(g1s), length(g2s), length(g3s)))
              surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
              surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
              pvalue <- 1- pchisq (surv.diff$chisq[1], df=2)
              plot(surv.fit, xlab="Time", ylab="Survival", main=paste("RNASeq2 -",myG), col=c(2,3,4))
              pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
              legTxt <- c("< 1st Q.", "Between 1st&3rd Q.", "> 3rd. Q.", pValueTxt)
              legend("topright", legend=legTxt, col=c(2,3,4,0), lty=c(1,1),cex=0.7)
            }
          }
        }
      }
      #}
    },
    "mRNAArray"=
    {
      for(jj in 1:length(dataObject@mRNAArray))
      {
        tmpMat1 <- dataObject@mRNAArray[[jj]]@DataMatrix
        sampleIDs1 <- colnames(tmpMat1)
        sampleIDs1 <- gsub(pattern="\\.",replacement="-",sampleIDs1)
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
        sampleIDs11 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
        colnames(tmpMat1) <- sampleIDs1
        tmpMat1 <- tmpMat1[,!duplicated(sampleIDs11)]
        colnames(tmpMat1) <- sampleIDs11[!duplicated(sampleIDs11)]
        for(myG in geneSymbols)
        {
          if(!is.na(any(match(rownames(tmpMat1),myG))) & any(match(rownames(tmpMat1),myG)))
          {
            tmpGeneExp <- as.numeric(tmpMat1[myG,])
            names(tmpGeneExp) <- colnames(tmpMat1)
            commonSamples <- intersect(rownames(sampleTimeCensor),names(tmpGeneExp))
            if(length(commonSamples) > 9)
            {
              if(numberofGroups==2)
              {
                g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[3]]
                g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[3]]
                time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),2])
                censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s),3])
                surv.group <- rep (1:2, c(length(g1s), length(g2s)))
                surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
                surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
                pvalue <- 1- pchisq (surv.diff$chisq[1], df=1)
                myFN <- gsub(x=dataObject@mRNAArray[[jj]]@Filename,pattern="__",replacement="@")
                myFN <- gsub(x=myFN,pattern="_",replacement="@")
                myFN <- paste(unlist(strsplit(x=myFN,split="@"))[3],unlist(strsplit(x=myFN,split="@"))[4],unlist(strsplit(x=myFN,split="@"))[5])
                plot(surv.fit, xlab="Time", ylab="Survival", main=paste(myFN,myG,sep="-"), col=c(2,4))
                pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
                legTxt <- c("< Median", ">= Median", pValueTxt)
                legend("topright", legend=legTxt, col=c(2,4,0), lty=c(1,1),cex=0.7)
              }
              else if(numberofGroups==3)
              {
                g1s <- names(tmpGeneExp)[tmpGeneExp < summary(tmpGeneExp)[2]]
                g2s <- names(tmpGeneExp)[tmpGeneExp >= summary(tmpGeneExp)[2] & tmpGeneExp <= summary(tmpGeneExp)[5]]
                g3s <- names(tmpGeneExp)[tmpGeneExp > summary(tmpGeneExp)[5]]
                time.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),2])
                censor.group <- as.numeric(sampleTimeCensor[c(g1s,g2s,g3s),3])
                surv.group <- rep (1:3, c(length(g1s), length(g2s), length(g3s)))
                surv.fit <- survfit (Surv(time.group,censor.group)~surv.group)
                surv.diff <- survdiff (Surv(time.group,censor.group)~surv.group)
                pvalue <- 1- pchisq (surv.diff$chisq[1], df=2)
                myFN <- gsub(x=dataObject@mRNAArray[[jj]]@Filename,pattern="__",replacement="@")
                myFN <- gsub(x=myFN,pattern="_",replacement="@")
                myFN <- paste(unlist(strsplit(x=myFN,split="@"))[3],unlist(strsplit(x=myFN,split="@"))[4],unlist(strsplit(x=myFN,split="@"))[5])
                plot(surv.fit, xlab="Time", ylab="Survival", main=paste(myFN,myG,sep="-"), col=c(2,3,4))
                pValueTxt <- paste ("p-value= ", format.pval(pvalue,digits=2), sep="")
                legTxt <- c("< 1st Q.", "Between 1st&3rd Q.", "> 3rd. Q.", pValueTxt)
                legend("topright", legend=legTxt, col=c(2,3,4,0), lty=c(1,1),cex=0.7)
              }
            }
          }
        }
      }
    })
  }
}
