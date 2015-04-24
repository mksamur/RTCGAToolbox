#' Extract data from \code{FirehoseData} object into S4 \code{ExpressionSet}
#' 
#' This function serves to extract and reorganize data into a structured S4 
#' object for analysis. An option is available to retreive clinical data 
#' and return an ExpressionSet class object (see \code{\link{ExpressionSet}})
#' otherwise, a data frame is returned. 
#' 
#' @param object A \code{FirehoseData} object from which to extract data. 
#' @param type The type of data to extract from the "FirehoseData" object.
#' @param phenoData Logical (default TRUE) includes additional clinic data, if available.
#' @return An \code{\link{ExpressionSet}} object for the selected data type or data frame if no phenoData is requested. 
#' Choices include: "RNAseq_Gene", "Clinic", "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA", "GISTIC_A", "GISTIC_T".
#' The "GISTIC_A" type of dataset represents GISTIC data by all genes. "GISTIC_T"" represents data thresholded by genes.
#' 
#' @examples 
#' 
#' \dontrun{
#' b2 <- extract(a2, "Methylation", phenoData=TRUE)
#' }
#' 
extract <- function(object, type, phenoData = TRUE){
  choices <- c("RNAseq_Gene", "miRNASeq_Gene",
               "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
               "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
               "miRNA_Array", "RPPA") 
  if(type %in% choices){
    slotreq <- grep(paste0("^", gsub("_", "", type))  , slotNames(object), 
                    ignore.case=TRUE, perl=TRUE, value=TRUE)
    elemlength <- length(getElement(object, slotreq))
    if(is(getElement(object, slotreq), "list")){
      if(elemlength > 1){
        dm <- lapply(getElement(object, slotreq), function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(dm, ncol))
        dm <- dm[[keeplist]]
        warning(paste("Taking the array platform with the greatest number of samples:", keeplist))
      } else if(elemlength == 1){
        dm <- getElement(object, slotreq)[[1]]@DataMatrix
      } else if(elemlength == 0){
        dm <- matrix(NA, nrow=0, ncol=0)
      }
    } else {
      dm <- getElement(object, slotreq)
    }
  } else if(type %in% c("GISTIC_A", "GISTIC_T")) {
    if(type=="GISTIC_A"){
      dm <- getElement(object@GISTIC, "AllByGene")
    } else {
      dm <- getElement(object@GISTIC, "ThresholedByGene")
    }
  } else {
    stop(paste("Type not yet supported or could not be matched."))
  }
  
  bcID <- function(barcodes, sample=FALSE, portion = FALSE, center=FALSE, collapse=FALSE){
    barcodes <- gsub("\\.", "-", barcodes)
    if(center){
      if(!sample){
        bcc <- sapply(strsplit(barcodes, "-"),"[", c(6:7))
        bcc <- tolower(t(bcc))
        return(bcc)
      } else {
        bcc <- sapply(strsplit(barcodes, "-"), "[", c(4:7))
        bcc <- tolower(t(bcc))
      }
    } else {
      bcc <- lapply(strsplit(barcodes, "-"), "[", c(1:3))
      bcc <- tolower(sapply(bcc, paste, collapse="-"))
      if(!sample & portion){
        port <- tolower(sapply(strsplit(barcodes, "-"), "[", 5))
        return(port)
      }
      if(sample & collapse){
        samp <- lapply(strsplit(barcodes, "-"), "[", c(1:4))
        samp <- tolower(sapply(samp, paste, collapse="-"))
        samp <- substr(samp, 1, (nchar(samp)-1))
        return(samp)
      }else if(sample & !collapse){
        samp <- tolower(sapply(strsplit(barcodes, "-"), "[", 4))
        samp <- as.numeric(substr(samp, 1, nchar(samp)-1))
        return(samp)
      }
    }
    return(bcc)
  }
  
  if(dim(dm)[1] == 0 | dim(dm)[2] == 0){
    stop("There is no data for that data type!")
  } else {
      if(slotreq=="Methylation"){
      annote <- dm[,c("Gene_Symbol", "Chromosome", "Genomic_Coordinate")]      
      dm <- apply(dm[4:ncol(dm)], 2, as.numeric, as.matrix)
      rownames(dm) <- rownames(dm)
      }
      dups <- colnames(dm)[duplicated(bcID(colnames(dm)))]
      
      # technical replicates and tumor/normals present
      if(length(dups) != 0){
        # from TCGA Code Tables Report  
        
        sample_type <- samptab[,2][match(bcID(colnames(dm), sample=TRUE), samptab[,1])]
        rightbc <- data.frame(sample_type, 
                              sample_code = substr(bcID(colnames(dm), sample=TRUE), 1,2), vial = substr(bcID(colnames(dm), sample=TRUE, center=TRUE)[,1], 3,3),
                              portion = substr(bcID(colnames(dm), portion=TRUE), 1,2), analyte = substr(bcID(colnames(dm), portion=TRUE), 3,3), 
                              plate = bcID(colnames(dm), center=T)[,1], center = bcID(colnames(dm), center=T)[,2],  stringsAsFactors=FALSE)
        rightbc[, "sample_code"] <- as.numeric(rightbc[,"sample_code"])
        rownames(rightbc) <- bcID(colnames(dm), sample=TRUE, collapse=TRUE)
              
        repeated <- dm[, bcID(colnames(dm)) %in% bcID(dups)]
        samps <- bcID(colnames(repeated), sample=TRUE)
        normals <- repeated[, samps >= 10 & samps <= 19]
        
        if(range(samps)[2] > 19){
          controls <- repeated[, samps >= 20 & samps <= 29]
          replicates <- repeated[, !bcID(colnames(repeated)) %in% bcID(colnames(normals)) & 
                                   !bcID(colnames(repeated)) %in% bcID(colnames(controls))]
        } else {      
        replicates <- repeated[, !bcID(colnames(repeated)) %in% bcID(colnames(normals))]
        }
        
        duplic <- colnames(replicates)[duplicated(bcID(colnames(replicates)))]
        d <- c()
        for (cc in seq(duplic)){
          d <- cbind(d, apply(replicates[,bcID(colnames(replicates)) %in% bcID(duplic[cc])], 1, mean))
        }
        colnames(d) <- bcID(duplic)
        dm <- cbind(dm[,!(bcID(colnames(dm)) %in% bcID(dups))], d)
        colnames(dm) <- bcID(colnames(dm))
      }
      
      if(phenoData){
        
        pd <- getElement(object, "Clinical")
        rownames(pd) <- gsub("\\.", "-", rownames(pd))
        if(length(pd)==0){
          stop("No clinical data available!")
        }
      
      npd <- pd[na.omit(match(bcID(colnames(dm)), rownames(pd))),]
      npd <- cbind(npd, rightbc[match(rownames(npd), bcID(rownames(rightbc))),])
      npd <-  npd[, -1]
      
      # getLinks(".Merge_Clinical.Level_1", "*.tar[.]gz$")
      cl_url <- "http://gdac.broadinstitute.org/runs/stddata__"
      cl_url <- paste0(cl_url,substr(getElement(object, "runDate"),1,4),"_",substr(getElement(object, "runDate"),5,6),"_",substr(getElement(object, "runDate"),7,8),"/data/")
      cl_url <- paste0(cl_url,getElement(object, "Dataset"),"/",getElement(object, "runDate"),"/")
      cl_url <- paste0(cl_url, "gdac.broadinstitute.org_", getElement(object, "Dataset"), ".Merge_Clinical.Level_1.", getElement(object, "runDate"), "00.0.0.tar.gz")
      download.file(url=cl_url, destfile=paste0(getElement(object, "Dataset"), "-ExClinical.tar.gz"), method="auto", quiet=TRUE, mode="w")
      fileList <- untar(paste0(getElement(object, "Dataset"), "-ExClinical.tar.gz"), list=TRUE)
      fileList <- grep(".clin.merged.txt", fileList, fixed = TRUE, value=TRUE)
      untar(paste0(getElement(object, "Dataset"),"-ExClinical.tar.gz"),files=fileList)
      file.rename(from=fileList,to=paste0(getElement(object, "runDate"),"-",getElement(object, "Dataset"),"-ExClinical.txt"))
      file.remove(paste0(getElement(object, "Dataset"),"-ExClinical.tar.gz"))
      unlink(strsplit(fileList[1],"/")[[1]][1], recursive = TRUE)
      extracl <- fread(paste0(getElement(object, "runDate"),"-",getElement(object, "Dataset"),"-ExClinical.txt"), data.table=FALSE)
      colnames(extracl)[-1] <- extracl[grep("patient_barcode", extracl[, 1]),][-1]
      rownames(extracl) <- extracl[, 1]      
      extracl <- extracl[,-1]
      extracl <- t(extracl)
      extracl <- extracl[,!grepl("patient_barcode", colnames(extracl))]
      
      phenoD <- merge(npd, extracl, "row.names")
      rownames(phenoD) <- phenoD[,"Row.names"]
      phenoD <- phenoD[, -1]
      
      ndm <- dm[,na.omit(match(rownames(phenoD), colnames(dm)))]
      
      
      if(identical(all.equal(rownames(phenoD), colnames(ndm)), TRUE)){
        eset <- ExpressionSet(ndm, AnnotatedDataFrame(phenoD))
        if(exists("annote")){
          featureData(eset) <- AnnotatedDataFrame(annote)
        }
        return(eset)
      }else{
        stop("Couldn't match up rownames of phenodata to colnames of numeric data")
        return(dm)
      }
  
      return(eset)
      
      } else {
        return(dm)
      }
      
    }
  }

  