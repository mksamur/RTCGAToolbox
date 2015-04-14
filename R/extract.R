#' Extract data from \code{FirehoseData} object into S4 \code{ExpressionSet}
#' 
#' This function serves to extract and reorganize data into a structured S4 
#' object for genomic analysis. An option is available to retreive additional
#' clinical data and include it as featureData. See \code{\link{ExpressionSet}}
#' 
#' @param object A \code{FirehoseData} object from which to extract data. 
#' @param type The type of data to extract from the "FirehoseData" object.
#' @param phenoData Logical (default TRUE) includes additional clinic data, if available.
#' @return An \code{\link{ExpressionSet}} object for the selected data type. 
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
        output <- lapply(getElement(object, slotreq), function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the array platform with the greatest number of samples:", keeplist))
      } else if(elemlength == 1){
        output <- getElement(object, slotreq)[[1]]@DataMatrix
      } else if(elemlength == 0){
        output <- matrix(NA, nrow=0, ncol=0)
      }
    } else {
      output <- getElement(object, slotreq)
    }
  } else if(type %in% c("GISTIC_A", "GISTIC_T")) {
    if(type=="GISTIC_A"){
      output <- getElement(object@GISTIC, "AllByGene")
    } else {
      output <- getElement(object@GISTIC, "ThresholedByGene")
    }
  } else {
    stop(paste("Type not yet supported or could not be matched."))
  }
  
  bcID <- function(barcodes, center=FALSE, sample=FALSE, collapse=FALSE){
    barcodes <- gsub("\\.", "-", barcodes)
    if(center){
      if(!sample){
        bcc <- sapply(strsplit(barcodes, "-"),"[", c(5:7))
        bcc <- tolower(t(bcc))
        return(bcc)
      } else {
        bcc <- sapply(strsplit(barcodes, "-"), "[", c(4:7))
        bcc <- tolower(t(bcc))
      }
    } else {
      bcc <- lapply(strsplit(barcodes, "-"), "[", c(1:3))
      bcc <- tolower(sapply(bcc, paste, collapse="-"))
      if(sample & collapse){
        samp <- lapply(strsplit(barcodes, "-"), "[", c(1:4))
        samp <- tolower(sapply(samp, paste, collapse="-"))
        samp <- substr(samp, 1, nchar(samp)[1]-1)
        return(samp)
      }else if(sample & !collapse){
        samp <- tolower(sapply(strsplit(barcodes, "-"), "[", 4))
        samp <- substr(samp, 1,2)
        return(samp)
      }
    }
    return(bcc)
  }
  
  if(dim(output)[1] == 0 | dim(output)[2] == 0){
    message("There is no data for that data type!")
  } else {
    if(phenoData){
            
      dm <- apply(output[4:ncol(output)], 2, as.numeric, as.matrix)
      rownames(dm) <- rownames(output)
      
      dups <- colnames(dm)[duplicated(bcID(colnames(dm)))]
      
      # technical replicates and tumor/normals present
      if(length(dups) != 0){
        repeated <- dm[, bcID(colnames(dm)) %in% bcID(dups)]
        samps <- as.numeric(bcID(colnames(repeated), sample=TRUE))
        normals <- repeated[, samps >= 10 & samps <= 19]
        controls <- repeated[, samps >= 20 & samps <= 29]
        # follow up on controls
        repeated <- repeated[, !bcID(colnames(repeated)) %in% bcID(colnames(normals))]
        duplic <- colnames(repeated)[duplicated(bcID(colnames(repeated)))]
        d <- c()
        for (cc in seq(duplic)){
          d <- cbind(d, apply(repeated[,bcID(colnames(repeated)) %in% bcID(duplic[cc])], 1, mean))
        }
        colnames(d) <- duplic
        dm <- cbind(dm[,!(bcID(colnames(dm)) %in% bcID(dups))], d)
      }
      
      # tumors <- dm[, samps <= 9]
      centerbc <- bcID(colnames(dm), center=TRUE, sample=TRUE)
      colnames(centerbc) <- c("sample", "portion", "plate", "center")
      rownames(centerbc) <- bcID(colnames(dm))
      
      pd <- getElement(object, "Clinical")
      rownames(pd) <- gsub("\\.", "-", rownames(pd))
      if(length(pd)==0){
        stop("No clinical data available!")
      }
      
      colnames(dm) <- bcID(colnames(dm))
      
      npd <- pd[na.omit(match(bcID(colnames(dm)), rownames(pd))),]
      npd <- merge(npd, centerbc, "row.names")
      rownames(npd) <- npd[, "Row.names"]
      npd <-  npd[, -1]
            
      ndm <- dm[,na.omit(match(rownames(npd), colnames(dm)))]
      
      if(identical(all.equal(rownames(npd), colnames(ndm)), TRUE)){
        eset <- ExpressionSet(ndm, AnnotatedDataFrame(npd))
      }else{
        stop("Couldn't match up rownames of phenodata to colnames of numeric data")
      }
      
      return(eset)
      #       ovextraclin <- fread("extraclinfilehere")
      #       ovextraclin[grep("patient_barcode", ovextraclin[, 1]),][-1]
      #       colnames(ovextraclin)[-1] <- ovextraclin[grep("patient_barcode", ovextraclin[, 1]),][-1]
      #       rownames(ovextraclin) <- ovextraclin[, 1]      
      #       jj <- ovextraclin[,match(rownames(npd), colnames(ovextraclin))]
      #       featureData(eset) <- AnnotatedDataFrame(jj)
      #       return(eset)
    } else {
      if(slotreq=="Methylation")
      {
        output <- output[, -c(1:3)]
        colnames(output) <- bcID(colnames(output), sample = TRUE, collapse=TRUE)
        return(output)  
      }
    }
  }
}
  