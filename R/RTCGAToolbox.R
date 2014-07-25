#Set classes
setClass("FirehoseCGHArray", representation(Filename = "character", DataMatrix = "data.frame"))
setClass("FirehoseMethylationArray", representation(Filename = "character", DataMatrix = "data.frame"))
setClass("FirehosemRNAArray", representation(Filename = "character", DataMatrix = "data.frame"))
setClass("FirehoseGISTIC", representation(Dataset = "character", AllByGene = "data.frame",ThresholedByGene="data.frame"))

setClass("FirehoseData", representation(Dataset = "character", Clinical = "data.frame", RNASeqGene = "matrix",
                                        RNASeq2GeneNorm="matrix",miRNASeqGene="matrix",CNASNP="data.frame",
                                        CNVSNP="data.frame",CNAseq="data.frame",CNACGH="list",Methylation="list",
                                        mRNAArray="list",miRNAArray="list",RPPAArray="list",Mutations="data.frame",
                                        GISTIC="FirehoseGISTIC"))

setClass("DGEResult", representation(Dataset = "character", Toptable = "data.frame"))
setClass("CorResult", representation(Dataset = "character", Correlations = "data.frame"))


#Main data client function
getFirehoseData <- function(dataset, runDate=NULL, gistic2_Date=NULL, RNAseq_Gene=FALSE,Clinic=TRUE,
                            miRNASeq_Gene=FALSE, RNAseq2_Gene_Norm=FALSE,
                            CNA_SNP=FALSE,CNV_SNP=FALSE,
                            CNA_Seq=FALSE,CNA_CGH=FALSE,Methylation=FALSE,Mutation=FALSE,mRNA_Array=FALSE,
                            miRNA_Array=FALSE,RPPA=FALSE,RNAseqNorm="raw_counts",RNAseq2Norm="normalized_count")
{
  #check parameters
  if(!class(dataset)=="character" || is.null(dataset) || !length(dataset) == 1 || nchar(dataset) < 2)
  {stop('Please set "dataset" parameter! You should specify one dataset name. Ex: dataset="BRCA"...')}
  runDatasets <- getFirehoseDatasets()
  if(!any(runDatasets==dataset)){stop('Please use valid dataset name! "getFirehoseDatasets" function gives you the vector of valid dataset names!')}
  
  
  if(!is.null(runDate))
  {
    if(!class(runDate)=="character" || !length(runDate) == 1 || !nchar(runDate) == 8)
    {stop('Please set "runDate" parameter! You should specify one Firehose run date. Ex: runDate="20140416"...')}
    
    runDateList <- getFirehoseRunningDates()
    if(!any(runDateList==runDate)){stop('Please use valid run date! "getFirehoseRunningDates" function gives you the vector of valid dates!')}
  }
  
  if(!is.null(gistic2_Date))
  {
    if(!class(gistic2_Date)=="character" || !length(gistic2_Date) == 1 || !nchar(gistic2_Date) == 8)
    {stop('Please set "gistic2_Date" parameter! You should specify one Firehose run date. Ex: gistic2_Date="20140115"...')}
    
    runGisticDate <- getFirehoseAnalyzeDates()
    if(!any(runGisticDate==gistic2_Date)){stop('Please use valid analyze date for GISTIC! "getFirehoseAnalyzeDates" function gives you the vector of valid dates!')}
  }
  
  if(is.null(gistic2_Date) & is.null(runDate)){stop("Please specify run date or/and gistic date!")}
  
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  resultClass <- new("FirehoseData", Dataset = dataset)
  
  tryCatch({
  
  if(!is.null(runDate))
  {
    ##build URL for getting file links
    fh_url <- "http://gdac.broadinstitute.org/runs/stddata__"
    fh_url <- paste(fh_url,substr(runDate,1,4),"_",substr(runDate,5,6),"_",substr(runDate,7,8),"/data/",sep="")
    fh_url <- paste(fh_url,dataset,"/",runDate,"/",sep="")
    doc = htmlTreeParse(fh_url, useInternalNodes = T)
    
    #Download clinical data
    if(Clinic)
    {
      keyWord = paste(dataset,".Clinical_Pick_Tier1.Level_4",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl("*.tar[.]gz$",plinks)]
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-Clinical.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-Clinical.tar.gz",sep=""),list=TRUE)
        fileList = fileList[grepl("*.clin.merged.picked.txt$",fileList)]
        untar(paste(dataset,"-Clinical.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-Clinical.txt",sep=""))
        file.remove(paste(dataset,"-Clinical.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        resultClass@Clinical <- read.delim(paste(dataset,"-Clinical.txt",sep=""),colClasses="character")
      }
    }
    
    #Download RNAseq gene level data
    if(RNAseq_Gene)
    {
      keyWord = paste("","Level_3__gene_expression__data.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl("*.Merge_rnaseq__.*._rnaseq__.*.tar[.]gz$",plinks)]
      for(i in trim(plinks))
      {
        
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-RNAseqGene.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-RNAseqGene.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]rnaseq__.*.__Level_3__gene_expression__data.data.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-RNAseqGene.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-RNAseqGene.txt",sep=""))
        file.remove(paste(dataset,"-RNAseqGene.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-RNAseqGene.txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 1:ncol(tmpCols)
        colOrder <- colOrder[tmpCols[1,] == RNAseqNorm]
        
        message("RNA-seq data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-RNAseqGene.txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        message(paste(nooflines,"genes data will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-RNAseqGene.txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"genes data has been imported!"))}
          else{message(paste((nooflines),"genes data has been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        tmpMat2 <- data.frame()
        for(i in 1:length(listMat))
        {
          if(i == 1){tmpMat2 = listMat[[i]]}
          else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
          listMat[[i]] <- 0 
        }
        tmpMat <- rbind(tmpMat2,tmpMat)
        rm(tmpMat2)
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- tmpMat[1,]
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        names1 <- tmpMat[,1]
        tmpMat <- tmpMat[,-1]
        tmpMat <- apply(tmpMat,2,as.numeric)
        rownames(tmpMat) <- names1
        
        resultClass@RNASeqGene <- tmpMat
        
      }
    }
    
    #Download RNAseq2 gene level data
    if(RNAseq2_Gene_Norm)
    {
      keyWord = paste("","Level_3__RSEM_genes_normalized__data.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl("*.Merge_rnaseqv2__.*._rnaseqv2__.*.tar[.]gz$",plinks)]
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]rnaseqv2__.*.__Level_3__RSEM_genes_normalized__data.data.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-RNAseq2GeneNorm.txt",sep=""))
        file.remove(paste(dataset,"-RNAseq2GeneNorm.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-RNAseq2GeneNorm.txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 1:ncol(tmpCols)
        colOrder <- colOrder[tmpCols[1,] == RNAseq2Norm]
        
        message("RNA-seq2 data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-RNAseq2GeneNorm.txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        message(paste(nooflines,"genes data will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-RNAseq2GeneNorm.txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"genes data has been imported!"))}
          else{message(paste((nooflines),"genes data has been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        tmpMat2 <- data.frame()
        for(i in 1:length(listMat))
        {
          if(i == 1){tmpMat2 = listMat[[i]]}
          else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
          listMat[[i]] <- 0 
        }
        tmpMat <- rbind(tmpMat2,tmpMat)
        rm(tmpMat2)
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- tmpMat[1,]
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        names1 <- tmpMat[,1]
        tmpMat <- tmpMat[,-1]
        tmpMat <- apply(tmpMat,2,as.numeric)
        rownames(tmpMat) <- names1
        
        resultClass@RNASeq2GeneNorm <- tmpMat
      }
    }
    
    #Download miRNAseq gene level data
    if(miRNASeq_Gene)
    {
      keyWord = paste("","Level_3__miR_gene_expression__data.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      message(plinks)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$",sep=""),plinks)]
      message(plinks)
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-miRNAseqGene.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-miRNAseqGene.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]mirnaseq__.*.__Level_3__miR_gene_expression__data.data.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-miRNAseqGene.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-miRNAseqGene.txt",sep=""))
        file.remove(paste(dataset,"-miRNAseqGene.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-miRNAseqGene.txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 1:ncol(tmpCols)
        colOrder <- colOrder[tmpCols[1,] == "read_count"]
        
        message("miRNA-seq data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-miRNAseqGene.txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        message(paste(nooflines,"genes data will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-miRNAseqGene.txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"genes data has been imported!"))}
          else{message(paste((nooflines),"genes data has been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        tmpMat2 <- data.frame()
        for(i in 1:length(listMat))
        {
          if(i == 1){tmpMat2 = listMat[[i]]}
          else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
          listMat[[i]] <- 0 
        }
        tmpMat <- rbind(tmpMat2,tmpMat)
        rm(tmpMat2)
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- tmpMat[1,]
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        names1 <- tmpMat[,1]
        tmpMat <- tmpMat[,-1]
        tmpMat <- apply(tmpMat,2,as.numeric)
        rownames(tmpMat) <- names1
        
        resultClass@miRNASeqGene <- tmpMat
      }
    }
    
    #Download CNA SNP data
    if(CNA_SNP)
    {
      keyWord = paste("","Level_3__segmented_scna_hg19__seg.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_snp__.*.__Level_3__segmented_scna_hg19__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-CNASNPHg19.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-CNASNPHg19.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]snp__.*.__Level_3__segmented_scna_hg19__seg.seg.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-CNASNPHg19.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-CNASNPHg19.txt",sep=""))
        file.remove(paste(dataset,"-CNASNPHg19.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpMat = read.delim(paste(dataset,"-CNASNPHg19.txt",sep=""),header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                            "numeric","numeric"))
        resultClass@CNASNP <- tmpMat
      }
    }
    
    #Download CNV SNP data
    if(CNV_SNP)
    {
      keyWord = paste("","Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_snp__.*.__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]snp__.*.__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-CNVSNPHg19.txt",sep=""))
        file.remove(paste(dataset,"-CNVSNPHg19.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpMat = read.delim(paste(dataset,"-CNVSNPHg19.txt",sep=""),header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                             "numeric","numeric"))
        resultClass@CNVSNP <- tmpMat
      }
    }
    
    #Download CNA DNAseq data
    if(CNA_Seq)
    {
      keyWord = paste("","__Level_3__segmentation__seg.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_cna__.*.dnaseq.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-CNAseq.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-CNAseq.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]cna__.*.__Level_3__segmentation__seg.seg.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-CNAseq.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-CNAseq.txt",sep=""))
        file.remove(paste(dataset,"-CNAseq.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpMat = read.delim(paste(dataset,"-CNAseq.txt",sep=""),header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                             "numeric","numeric"))
        resultClass@CNAseq <- tmpMat
      }
    }
    
    #Download CNA CGH data
    if(CNA_CGH)
    {
      keyWord = paste("","__Level_3__segmentation__seg.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_cna__.*.cgh.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$",sep=""),plinks)]
      
      dataLists <- list()
      listCount = 1
      for(i in trim(plinks))
      {
        download_link = paste(fh_url,i,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-CNACGH.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-CNACGH.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]cna__.*.__Level_3__segmentation__seg.seg.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-CNACGH.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-CNACGH-",listCount,".txt",sep=""))
        file.remove(paste(dataset,"-CNACGH.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpMat = read.delim(paste(dataset,"-CNACGH-",listCount,".txt",sep=""),header=TRUE,colClasses=c("character","numeric","numeric",
                                                                                         "numeric","numeric"))
        tmpReturn <- new("FirehoseCGHArray",Filename=i,DataMatrix=tmpMat)
        dataLists[[listCount]] <- tmpReturn
        listCount = listCount + 1
      }
      resultClass@CNACGH <- dataLists
    }
    
    #Download methylation
    if(Methylation)
    {
      keyWord = paste("","__Level_3__within_bioassay_data_set_function__data.Level_3",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks = xpathSApply(doc, keyWord, xmlValue)
      plinks = plinks[grepl(paste("*.",dataset,"[.]Merge_methylation__.*.methylation.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$",sep=""),plinks)]
      
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        download_link = paste(fh_url,ii,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-Methylation.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-Methylation.tar.gz",sep=""),list=TRUE)
        grepSearch = paste("*.",dataset,"[.]methylation__.*.__Level_3__within_bioassay_data_set_function__data.data.txt$",sep="")
        fileList = fileList[grepl(grepSearch,fileList)]
        untar(paste(dataset,"-Methylation.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-Methylation-",listCount,".txt",sep=""))
        file.remove(paste(dataset,"-Methylation.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-Methylation-",listCount,".txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 1:ncol(tmpCols)
        colOrder <- colOrder[tmpCols[1,] == "Beta_value"]
        
        message("Methylation data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-Methylation-",listCount,".txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        if(nooflines > 50000)
        {
          message(paste(ii,"won't be imported due to high data volume!"))
          break
        }
        message(paste(nooflines,"rows will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-Methylation-",listCount,".txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,3,4,5,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,3,4,5,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"rows have been imported!"))}
          else{message(paste((nooflines),"rows have been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        tmpMat2 <- data.frame()
        for(i in 1:length(listMat))
        {
          if(i == 1){tmpMat2 = listMat[[i]]}
          else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
          listMat[[i]] <- 0 
        }
        tmpMat <- rbind(tmpMat2,tmpMat)
        rm(tmpMat2)
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- c("CompositeElementREF","Gene_Symbol","Chromosome","Genomic_Coordinate",tmpMat[1,5:ncol(tmpMat)])
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        names1 <- tmpMat[,1]
        tmpMat <- tmpMat[,-1]
        #tmpMat <- apply(tmpMat,2,as.numeric)
        rownames(tmpMat) <- names1
        
        
        tmpReturn <- new("FirehoseMethylationArray",Filename=ii,DataMatrix=tmpMat)
        dataLists[[listCount]] <- tmpReturn
        listCount = listCount + 1
      }
      resultClass@Methylation <- dataLists
    }
    
    #Download mRNA array
    if(mRNA_Array)
    {
      keyWord = paste("","Merge_transcriptome__agilentg4502a_07",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks1 = xpathSApply(doc, keyWord, xmlValue)
      plinks1 = plinks1[grepl(paste("*.",dataset,"[.]Merge_transcriptome__agilentg4502a_.*.__Level_3__unc_lowess_normalization_gene_level__data.Level_3.*.tar[.]gz$",sep=""),plinks1)]
      
      keyWord = paste("","Merge_transcriptome__ht_hg_u133a",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks2 = xpathSApply(doc, keyWord, xmlValue)
      plinks2 = plinks2[grepl(paste("*.",dataset,"[.]Merge_transcriptome__ht_hg_u133a__.*.__Level_3__gene_rma__data.Level_3.*.tar[.]gz$",sep=""),plinks2)]
      
      keyWord = paste("","Merge_exon__huex_1_0_st_v2",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks3 = xpathSApply(doc, keyWord, xmlValue)
      plinks3 = plinks3[grepl(paste("*.",dataset,"[.]Merge_exon__huex_1_0_st_v2__.*.__Level_3__quantile_normalization_gene__data.Level_3.*.tar[.]gz$",sep=""),plinks3)]
      
      plinks = c(plinks1,plinks2,plinks3)
      plinks = unique(plinks[plinks != ""])
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        download_link = paste(fh_url,ii,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-mRNAArray.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-mRNAArray.tar.gz",sep=""),list=TRUE)
        grepSearch = "MANIFEST.txt"
        fileList = fileList[!grepl(grepSearch,fileList)]
        untar(paste(dataset,"-mRNAArray.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-mRNAArray-",listCount,".txt",sep=""))
        file.remove(paste(dataset,"-mRNAArray.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-mRNAArray-",listCount,".txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 2:ncol(tmpCols)
        #colOrder <- colOrder[tmpCols[1,] == "Signal"]
        
        message("mRNA data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-mRNAArray-",listCount,".txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        if(nooflines > 50000)
        {
          message(paste(ii,"won't be imported due to high data volume!"))
          break
        }
        message(paste(nooflines,"rows will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-mRNAArray-",listCount,".txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"rows have been imported!"))}
          else{message(paste((nooflines),"rows have been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        tmpMat2 <- data.frame()
        for(i in 1:length(listMat))
        {
          if(i == 1){tmpMat2 = listMat[[i]]}
          else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
          listMat[[i]] <- 0 
        }
        tmpMat <- rbind(tmpMat2,tmpMat)
        rm(tmpMat2)
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- c("Gene_Symbol",tmpMat[1,2:ncol(tmpMat)])
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        #names1 <- tmpMat[,1]
        #tmpMat <- tmpMat[,-1]
        #tmpMat2 <- apply(tmpMat[,2:ncol(tmpMat)],2,as.numeric)
        #tmpMat <- cbind(tmpMat[,1],tmpMat2)
        #colnames(tmpMat) <- c("Gene_Symbol",colnames(tmpMat)[2:ncol(tmpMat)])
        #rownames(tmpMat) <- names1
        
        
        tmpReturn <- new("FirehosemRNAArray",Filename=ii,DataMatrix=data.frame(tmpMat))
        dataLists[[listCount]] <- tmpReturn
        listCount = listCount + 1
      }
      resultClass@mRNAArray <- dataLists
    }
    
    
    #Download miRNA array
    if(miRNA_Array)
    {
      keyWord = paste("","h_mirna_8x15k",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks1 = xpathSApply(doc, keyWord, xmlValue)
      plinks1 = plinks1[grepl(paste("*.",dataset,"[.]Merge_mirna__h_mirna_8x15k.*.data.Level_3.*.tar[.]gz$",sep=""),plinks1)]
      
      plinks = plinks1#c(plinks1,plinks2,plinks3)
      plinks = unique(plinks[plinks != ""])
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        download_link = paste(fh_url,ii,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-miRNAArray.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-miRNAArray.tar.gz",sep=""),list=TRUE)
        grepSearch = "MANIFEST.txt"
        fileList = fileList[!grepl(grepSearch,fileList)]
        untar(paste(dataset,"-miRNAArray.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-miRNAArray-",listCount,".txt",sep=""))
        file.remove(paste(dataset,"-miRNAArray.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-miRNAArray-",listCount,".txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 2:ncol(tmpCols)
        #colOrder <- colOrder[tmpCols[1,] == "Signal"]
        
        message("mRNA data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-miRNAArray-",listCount,".txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        if(nooflines > 50000)
        {
          message(paste(ii,"won't be imported due to high data volume!"))
          break
        }
        message(paste(nooflines,"rows will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-miRNAArray-",listCount,".txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"rows have been imported!"))}
          else{message(paste((nooflines),"rows have been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        if(nooflines > 1000)
        {
          tmpMat2 <- data.frame()
          for(i in 1:length(listMat))
          {
            if(i == 1){tmpMat2 = listMat[[i]]}
            else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
            listMat[[i]] <- 0 
          }
          tmpMat <- rbind(tmpMat2,tmpMat)
          rm(tmpMat2)
        }
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- c("miRGene_Symbol",tmpMat[1,2:ncol(tmpMat)])
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        #names1 <- tmpMat[,1]
        #tmpMat <- tmpMat[,-1]
        #tmpMat2 <- apply(tmpMat[,2:ncol(tmpMat)],2,as.numeric)
        #tmpMat <- cbind(tmpMat[,1],tmpMat2)
        #colnames(tmpMat) <- c("Gene_Symbol",colnames(tmpMat)[2:ncol(tmpMat)])
        #rownames(tmpMat) <- names1
        
        
        tmpReturn <- new("FirehosemRNAArray",Filename=ii,DataMatrix=data.frame(tmpMat))
        dataLists[[listCount]] <- tmpReturn
        listCount = listCount + 1
      }
      resultClass@miRNAArray <- dataLists
    }
    
    #Download RPPA array
    if(RPPA)
    {
      keyWord = paste("","rppa_core",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks1 = xpathSApply(doc, keyWord, xmlValue)
      plinks1 = plinks1[grepl(paste("*.",dataset,"[.]Merge_protein_exp.*.protein_normalization__data.Level_3.*.tar[.]gz$",sep=""),plinks1)]
      
      plinks = plinks1#c(plinks1,plinks2,plinks3)
      plinks = unique(plinks[plinks != ""])
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        download_link = paste(fh_url,ii,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-RPPAArray.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-RPPAArray.tar.gz",sep=""),list=TRUE)
        grepSearch = "MANIFEST.txt"
        fileList = fileList[!grepl(grepSearch,fileList)]
        untar(paste(dataset,"-RPPAArray.tar.gz",sep=""),files=fileList)
        
        file.rename(from=fileList,to=paste(dataset,"-RPPAArray-",listCount,".txt",sep=""))
        file.remove(paste(dataset,"-RPPAArray.tar.gz",sep=""))
        delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
        
        #Get selected type only
        tmpCols = read.delim(paste(dataset,"-RPPAArray-",listCount,".txt",sep=""),nrows=1,colClasses="character")
        colOrder <- 2:ncol(tmpCols)
        #colOrder <- colOrder[tmpCols[1,] == "Signal"]
        
        message("mRNA data will be imported! This may take some times!")
        testcon <- file(paste(dataset,"-RPPAArray-",listCount,".txt",sep=""),open="r")
        readsizeof <- 1000
        nooflines <- 0
        ( while((linesread <- length(readLines(testcon,readsizeof))) > 0 )
          nooflines <- nooflines+linesread )
        close(testcon)
        message(nooflines)
        
        if(nooflines > 50000)
        {
          message(paste(ii,"won't be imported due to high data volume!"))
          break
        }
        message(paste(nooflines,"rows will be imported!"))
        
        tmpMat <- data.frame()
        inputfile<-file(paste(dataset,"-RPPAArray-",listCount,".txt",sep=""),open="r")
        listMat <- list()
        itemcount = 1
        
        N = nooflines
        chunksize<-100
        nchunks<- ceiling(N/100)
        
        for(i in 1:nchunks){
          chunk<-read.delim(inputfile,nrows=chunksize,colClasses="character",header=FALSE)
          if(i == 1)
          {
            tmpMat <- chunk[,c(1,colOrder)]
          }
          else
          {
            tmpMat <- rbind(tmpMat,chunk[,c(1,colOrder)])
          }
          if((chunksize*i) < nooflines){message(paste((chunksize*i),"rows have been imported!"))}
          else{message(paste((nooflines),"rows have been imported!"))}
          
          if( ((chunksize*i) %% 1000) ==0 )
          {
            listMat[[itemcount]] <- tmpMat
            tmpMat <- data.frame()
            itemcount = itemcount + 1
          }
        }
        
        if(nooflines > 1000)
        {
          tmpMat2 <- data.frame()
          for(i in 1:length(listMat))
          {
            if(i == 1){tmpMat2 = listMat[[i]]}
            else{tmpMat2 = rbind(tmpMat2,listMat[[i]])}
            listMat[[i]] <- 0 
          }
          tmpMat <- rbind(tmpMat2,tmpMat)
          rm(tmpMat2)
        }
        
        #close(inputfile)
        closeAllConnections()
        colnames(tmpMat) <- c("Protein_Symbol",tmpMat[1,2:ncol(tmpMat)])
        tmpMat <- tmpMat[-c(1:2),]
        removeQM <- grepl("\\?\\|",tmpMat[,1])
        tmpMat <- tmpMat[!removeQM,]
        #names1 <- tmpMat[,1]
        #names2 <- sapply(names1,function(x){unlist(strsplit(x,"\\|"))[1]})
        #names1 <- tmpMat[,1]
        #tmpMat <- tmpMat[,-1]
        #tmpMat2 <- apply(tmpMat[,2:ncol(tmpMat)],2,as.numeric)
        #tmpMat <- cbind(tmpMat[,1],tmpMat2)
        #colnames(tmpMat) <- c("Gene_Symbol",colnames(tmpMat)[2:ncol(tmpMat)])
        #rownames(tmpMat) <- names1
        
        
        tmpReturn <- new("FirehosemRNAArray",Filename=ii,DataMatrix=data.frame(tmpMat))
        dataLists[[listCount]] <- tmpReturn
        listCount = listCount + 1
      }
      resultClass@RPPAArray <- dataLists
    }
    
    
    #Download RPPA array
    if(Mutation)
    {
      keyWord = paste("","Mutation_Packager_Calls",sep="")
      keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
      plinks1 = xpathSApply(doc, keyWord, xmlValue)
      plinks1 = plinks1[grepl(paste("*.",dataset,"[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$",sep=""),plinks1)]
      
      plinks = plinks1#c(plinks1,plinks2,plinks3)
      plinks = unique(plinks[plinks != ""])
      dataLists <- list()
      listCount = 1
      for(ii in trim(plinks))
      {
        download_link = paste(fh_url,ii,sep="")
        download.file(url=download_link,destfile=paste(dataset,"-Mutation.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
        fileList <- untar(paste(dataset,"-Mutation.tar.gz",sep=""),list=TRUE)
        grepSearch = "MANIFEST.txt"
        fileList = fileList[!grepl(grepSearch,fileList)]
        
        ###
        untar(paste(dataset,"-Mutation.tar.gz",sep=""),files=fileList)
        retMutations <- do.call("rbind",lapply(fileList,FUN=function(files){
          read.delim(files,header=TRUE,colClasses="character")
        }))
        delFodler <- paste(getwd(),"/",strsplit(fileList[1],"/")[[1]][1],sep="")
        unlink(delFodler, recursive = TRUE)
        file.remove(paste(dataset,"-Mutation.tar.gz",sep=""))
        ###
        
        #myMutFiles <- list()
        #countPos=1
        #for(myFiles in fileList)
        #{
        #  untar(paste(dataset,"-Mutation.tar.gz",sep=""),files=myFiles)
        #  tmpCols = read.delim(myFiles,header=TRUE,colClasses="character")
        #  myMutFiles[[countPos]] <- tmpCols
        #  countPos = countPos + 1
        #  
        #  delFodler <- paste(getwd(),"/",strsplit(myFiles,"/")[[1]][1],sep="")
        #  unlink(delFodler, recursive = TRUE)
        #}
        #file.remove(paste(dataset,"-Mutation.tar.gz",sep=""))
        #
        #retMutations <- data.frame()
        #for(i in 1:length(myMutFiles))
        #{
        #  if(i == 1){retMutations = myMutFiles[[i]]}
        #  else
        #  {
        #    retMutations = rbind(retMutations,myMutFiles[[i]])
        #  }
        #  
        #}
        
        write.table(retMutations,file=paste(dataset,"-Mutations-AllSamples.txt",sep=""),sep="\t",row.names=F,quote=F)
        resultClass@Mutations <- retMutations
      }
      
    }
    
    
    
    
  }
  
  if(!is.null(gistic2_Date))
  {
    ##build URL for getting file links
    fh_url <- "http://gdac.broadinstitute.org/runs/analyses__"
    fh_url <- paste(fh_url,substr(gistic2_Date,1,4),"_",substr(gistic2_Date,5,6),"_",substr(gistic2_Date,7,8),"/data/",sep="")
    fh_url <- paste(fh_url,dataset,"/",gistic2_Date,"/",sep="")
    doc = htmlTreeParse(fh_url, useInternalNodes = T)
    
    keyWord = paste("","CopyNumber_Gistic2.Level_4",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = xpathSApply(doc, keyWord, xmlValue)
    plinks = plinks[grepl(paste("*.",dataset,"-TP[.]CopyNumber_Gistic2[.]Level_4.*.tar[.]gz$",sep=""),plinks)]
    
    for(ii in trim(plinks))
    {
      download_link = paste(fh_url,ii,sep="")
      download.file(url=download_link,destfile=paste(dataset,"-Gistic2.tar.gz",sep=""),method="auto",quiet = FALSE, mode = "w")
      
      fileList <- untar(paste(dataset,"-Gistic2.tar.gz",sep=""),list=TRUE)
      grepSearch = "all_data_by_genes.txt"
      fileList = fileList[grepl(grepSearch,fileList)]
      untar(paste(dataset,"-Gistic2.tar.gz",sep=""),files=fileList)
      tmpCNAll = read.delim(fileList,header=TRUE,colClasses="character")
      file.rename(from=fileList,to=paste(dataset,"-all_data_by_genes.txt",sep=""))
      
      fileList <- untar(paste(dataset,"-Gistic2.tar.gz",sep=""),list=TRUE)
      grepSearch = "all_thresholded.by_genes.txt"
      fileList = fileList[grepl(grepSearch,fileList)]
      untar(paste(dataset,"-Gistic2.tar.gz",sep=""),files=fileList)
      tmpCNThreshhold = read.delim(fileList,header=TRUE,colClasses="character")
      file.rename(from=fileList,to=paste(dataset,"-all_thresholded.by_genes.txt",sep=""))
      
      
      delFodler <- paste(getwd(),"/",strsplit(fileList,"/")[[1]][1],sep="")
      unlink(delFodler, recursive = TRUE)
      
      
      tmpReturn <- new("FirehoseGISTIC",Dataset=dataset,AllByGene=data.frame(tmpCNAll),
                       ThresholedByGene=data.frame(tmpCNThreshhold))
      resultClass@GISTIC <- tmpReturn
      
    }
  }
  
  
  return(resultClass)
  
  
  },
  error=function(cond) {
    message("")
    message("")
    message("Dataset download error! Please check your connection and try again!")
    message("")
    message("")
    message("Detail Error Message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    message("")
    message("")
    message("Dataset download error! Please check your connection and try again!")
    message("If you keep getting this error message and you believe that your connection is stable you may use our Q&A Forum to get help from our community!")
    message("")
    message("")
    #message(cond)
    # Choose a return value in case of warning
    return(NULL)
  },
  finally={
    
  }) 
  
}
#####

#Check web resources
getFirehoseRunningDates <- function(last=NULL){
  check.integer <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
  }
  if(is.null(last)){
    runDate <- read.table("http://www.canevolve.org/fmineRdate.txt",
                          header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    return(runDate)
  }
  else if(check.integer(last))
  {
    runDate <- read.table("http://www.canevolve.org/fmineRdate.txt",
                             header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    if(last < length(runDate)){runDate <- runDate[1:last]}
    return(runDate)
  }
  else
  {
    stop('"last" must be integer')
  }
}

getFirehoseDatasets <- function(){
  runDataset <- read.table("http://www.canevolve.org/fmineRdataset.txt",
                           header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  runDataset <- as.character(runDataset[,1])
  return(runDataset)
}

getFirehoseAnalyzeDates <- function(last=NULL){
  check.integer <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
  }
  if(is.null(last)){
    runDate <- read.table("http://www.canevolve.org/fmineRgistic.txt",
                          header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    return(runDate)
  }
  else if(check.integer(last))
  {
    runDate <- read.table("http://www.canevolve.org/fmineRgistic.txt",
                          header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    if(last < length(runDate)){runDate <- runDate[1:last]}
    return(runDate)
  }
  else
  {
    stop('"last" must be integer')
  }
}

#####

#Analyze functions
getDiffExpressedGenes <- function(dataObject,DrawPlots=TRUE,adj.method="BH",adj.pval=0.05,raw.pval=0.05,logFC=2,
                                  hmTopUpN=100,hmTopDownN=100)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  
  validMatrix <- character()
  #check expression data matrices
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
  if(length(dataObject@mRNAArray) > 0){validMatrix <- append(validMatrix,"mRNAArray")}
  
  if(length(validMatrix) == 0){stop("There is no valid expression data in the object!")}
  
  if(class(DrawPlots) != "logical" | is.null(DrawPlots)){stop("DrawPlots must be logical!")}
  
  if(is.null(adj.method) | is.na(adj.method) | (tmp %in% c("BH","BY","holm","none"))){adj.method="BH"}
  if(is.null(adj.pval) | is.na(adj.pval) | length(adj.pval) > 1 | adj.pval > 1 | adj.pval < 0){adj.pval=0.05}
  if(is.null(raw.pval) | is.na(raw.pval) | length(raw.pval) > 1 | raw.pval > 1 | raw.pval < 0){raw.pval=0.05}
  if(is.null(logFC) | is.na(logFC) | length(logFC) > 1 | logFC < 0 ){logFC=2}
  if(is.null(hmTopUpN) | is.na(hmTopUpN) | length(hmTopUpN) > 1 | hmTopUpN < 0){hmTopUpN=100}
  if(is.null(hmTopDownN) | is.na(hmTopDownN) | length(hmTopDownN) > 1 | hmTopDownN < 0){hmTopDownN=100}
  
  listResults <- list()
  
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  for(i in validMatrix)
  {
    if(i == "RNASeq")
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
          voomMat <- dataObject@RNASeqGene[meanCounts > 10,c(normalSamples,tumorSamples)]
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          v <- voom(voomMat,design,plot=DrawPlots)
          fit <- lmFit(v,design)
          cont.matrix <- makeContrasts(TumorvsNormal=Tumor-Normal, levels=design)
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
    }
    
    if(i == "RNASeq2")
    {
      chkTmp <- as.numeric(dataObject@RNAseq2_Gene_Norm[1,])
      if(all(is.wholenumber(chkTmp)) == FALSE){warning("RNASeq2 data does not look like raw counts! We will skip this data!")}
      else
      {
        sampleIDs <- colnames(dataObject@RNAseq2_Gene_Norm)
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
          meanCounts <- apply(dataObject@RNASeqGene,1,mean)
          voomMat <- dataObject@RNASeqGene[meanCounts > 10,c(normalSamples,tumorSamples)]
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          v <- voom(voomMat,design,plot=DrawPlots)
          fit <- lmFit(v,design)
          cont.matrix <- makeContrasts(TumorvsNormal=Tumor-Normal, levels=design)
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
    }
    
    if(i == "mRNAArray")
    {
      for(j in 1:length(dataObject@mRNAArray))
      {
        tmpObj <- dataObject@mRNAArray[[j]]
        genes <- tmpObj@DataMatrix[,1]
        genes <- setdiff(genes,genes[duplicated(genes)])
        geneMat <- tmpObj@DataMatrix[tmpObj@DataMatrix[,1] %in% genes,]
        rownames(geneMat) <- geneMat[,1]
        geneMat <- geneMat[,2:ncol(geneMat)]
        
        sampleIDs <- colnames(geneMat)[]
        samplesDat <- data.frame(matrix(nrow=length(sampleIDs),ncol=7))
        rownames(samplesDat) <- sampleIDs
        for(j in 1:length(sampleIDs))
        {
          if(grepl(".",sampleIDs[j]))
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
          suppressWarnings(geneMat <- apply(geneMat,2,as.numeric))
          rownames(geneMat) <- rN
          colnames(geneMat) <- cN
          design <- model.matrix (~0 + factor(c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))))
          colnames (design) <-c ("Normal", "Tumor")
          fit <- lmFit(geneMat,design)
          cont.matrix <- makeContrasts(TumorvsNormal=Tumor-Normal, levels=design)
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
                try(heatmap(v,col=bluered,scale="row",main=tmpObj@Filename,Colv=NA),silent=FALSE)
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
    }
    
    
  }
  
  
  
  return(listResults)
  
}

getCNGECorrelation <- function(dataObject,adj.method="BH",adj.pval=0.05,raw.pval=0.05)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  
  validMatrix <- character()
  #check expression data matrices
  if(length(dataObject@mRNAArray) > 0){validMatrix <- append(validMatrix,"mRNAArray")}
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
  
  if(is.null(adj.method) | is.na(adj.method) | (adj.method %in% c("BH","BY","holm","none"))){adj.method="BH"}
  if(is.null(adj.pval) | is.na(adj.pval) | length(adj.pval) > 1 | adj.pval > 1 | adj.pval < 0){adj.pval=0.05}
  if(is.null(raw.pval) | is.na(raw.pval) | length(raw.pval) > 1 | raw.pval > 1 | raw.pval < 0){raw.pval=0.05}
  
  
  if(length(validMatrix) == 0){stop("There is no valid expression data in the object!")}
  
  if(dim(dataObject@GISTIC@AllByGene)[1] == 0 | dim(dataObject@GISTIC@AllByGene)[2] == 0 ){stop("There is no GISTIC data!")}
  
  listResults <- list()
  
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
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
          rnaseqGenes <- rownames(tmpMat1)
          rnaseqGenes2 <- character()
          for(rg in rnaseqGenes)
          {
            rnaseqGenes2 <- append(rnaseqGenes2,as.character(strsplit(rg,"\\|")[[1]][1]))
          }
          
          rnaFrame <- data.frame(rnaseqGenes,rnaseqGenes2)
          rnaFrame <- rnaFrame[!duplicated(rnaFrame[,2]),]
          rnaseqGenes2 <- rnaFrame[,2]
          names(rnaseqGenes2) <- rnaFrame[,1]
          
          tmpMat1 <- tmpMat1[names(rnaseqGenes2),]
          rownames(tmpMat1) <- rnaseqGenes2
          
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
        
        sampleIDs1 <- sampleIDs1[2:length(sampleIDs1)]
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
        commonGenes <- intersect(dataObject@mRNAArray[[jj]]@DataMatrix[,1],dataObject@GISTIC@AllByGene[,1])
        
        if(length(commonSamples) > 5)
        {
          tmpMat1 <- dataObject@mRNAArray[[jj]]@DataMatrix
          rownames(tmpMat1) <- tmpMat1[,1]
          tmpMat1 <- tmpMat1[,2:ncol(tmpMat1)]
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

getSurvival <- function(dataObject,numberofGroups=2,geneSymbols,sampleTimeCensor)
{
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  
  validMatrix <- character()
  #check expression data matrices
  if(dim(dataObject@RNASeqGene)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq")}
  if(dim(dataObject@RNASeq2GeneNorm)[1] > 0 & dim(dataObject@RNASeqGene)[2] > 0){validMatrix <- append(validMatrix,"RNASeq2")}
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
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  for(i in validMatrix)
  {    
    if(i == "RNASeq")
    {
      chkTmp <- as.numeric(dataObject@RNASeqGene[1,])
      controlVal = FALSE
      if(all(is.wholenumber(chkTmp)) == TRUE)
      {
        #warning("Current version of survival tool only works with normalized RNASeq data!")
        controlVal=TRUE
      }
      #else
      #{
        
        tmpMat1 <- dataObject@RNASeqGene
        rnaseqGenes <- rownames(tmpMat1)
        rnaseqGenes2 <- character()
        for(rg in rnaseqGenes)
        {
          rnaseqGenes2 <- append(rnaseqGenes2,as.character(strsplit(rg,"\\|")[[1]][1]))
        }
        
        rnaFrame <- data.frame(rnaseqGenes,rnaseqGenes2)
        rnaFrame <- rnaFrame[!duplicated(rnaFrame[,2]),]
        rnaseqGenes2 <- rnaFrame[,2]
        names(rnaseqGenes2) <- rnaFrame[,1]
        
        tmpMat1 <- tmpMat1[names(rnaseqGenes2),]
        rownames(tmpMat1) <- rnaseqGenes2
        
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
    }
    if(i == "mRNAArray")
    {
      for(jj in 1:length(dataObject@mRNAArray))
      {
        tmpMat1 <- dataObject@mRNAArray[[jj]]@DataMatrix
        rownames(tmpMat1) <- tmpMat1[,1]
        tmpMat1 <- tmpMat1[,2:ncol(tmpMat1)]
        
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
    }
  }
  
}

getReport <- function(dataObject,DGEResult1=NULL,DGEResult2=NULL,geneLocations)
{

  
  if(is.null(dataObject) | class(dataObject) != "FirehoseData")
  {stop("Please set a valid object! dataObject must be set as FirehoseData class!")}
  
  
  if(!is.null(DGEResult1) & class(DGEResult1) != "DGEResult"){stop("DGEResult1 must be DGEResult class!")}
  if(!is.null(DGEResult2) & class(DGEResult2) != "DGEResult"){stop("DGEResult2 must be DGEResult class!")}

  
  pdf(file=paste(dataObject@Dataset,"-reportImage.pdf",sep=""),height=30,width=30)
  plotpos = 1;
  require("RCircos")
  data(UCSC.HG19.Human.CytoBandIdeogram)
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
  
  if(!is.null(dataObject@GISTIC) & class(dataObject@GISTIC)=="FirehoseGISTIC")
  {
    cnMat <- dataObject@GISTIC@ThresholedByGene
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
      message(length(cnMat2[cnMat2 > 0.3 | cnMat2 < -0.3]))
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
    message(length(intGenes))
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
#####