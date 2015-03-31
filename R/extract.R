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
  typematch <- match.arg(type,
                         choices=c("RNAseq_Gene", "miRNASeq_Gene",
                                   "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
                                   "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
                                   "miRNA_Array", "RPPA", "GISTIC_A", "GISTIC_T"))
  if(identical(typematch, "RNAseq_Gene")){
    output <- getElement(object, "RNASeqGene")
  }else if(identical(typematch, "RNAseq2_Gene_Norm")){
    output <- getElement(object, "RNASeq2GeneNorm")
  }else if(identical(typematch, "miRNASeq_Gene")){
    output <- getElement(object, "miRNASeqGene")
  }else if(identical(typematch, "CNA_SNP")){
    output <- getElement(object, "CNASNP")
  }else if(identical(typematch, "CNV_SNP")){
    output <- getElement(object, "CNVSNP")
  }else if(identical(typematch, "CNA_Seq")){
    output <- getElement(object, "CNAseq")
  }else if(identical(typematch, "CNA_CGH")){
    if(is(object@CNACGH, "list")){
      if(length(object@CNACGH) > 1){
        output <- lapply(object@CNACGH, function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the CNACGH array platform with the greatest number of samples:", keeplist))
      }else if(length(object@CNACGH) == 1){
        output <- object@CNACGH[[1]]@DataMatrix
      }else{
        output <- matrix(NA, nrow=0, ncol=0)
      }
    }
  }else if(identical(typematch, "Mutation")){
    output <- getElement(object, "Mutations")
  }else if(identical(typematch, "RPPA")){
    if(is(object@RPPAArray, "list")){
      if(length(object@RPPAArray) > 1){
        output <- lapply(object@RPPAArray, function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the RPPA array platform with the greatest number of samples:", keeplist))
      }else if(length(object@RPPAArray) == 1){
        output <- object@RPPAArray[[1]]@DataMatrix
      }else{
        output <- matrix(NA, nrow=0, ncol=0)
      }
    }
  }else if(identical(typematch, "Methylation")){
    if(is(object@Methylation, "list")){
      if(length(object@Methylation) > 1){
        output <- lapply(object@Methylation, function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the Methylation array platform with the greatest number of samples:", keeplist))
      }else if(length(object@Methylation)==1){
        output <- object@Methylation[[1]]@DataMatrix
      }else{
        output <- matrix(NA, nrow=0, ncol=0)
      }
    }
  }else if(identical(typematch, "miRNA_Array")){
    if(is(object@miRNAArray, "list")){
      if(length(object@miRNAArray)>1){
        output <- lapply(object@miRNAArray, function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the miRNA array platform with the greatest number of samples:", keeplist))
      }else if(length(object@miRNAArray)==1){
        output <- object@miRNAArray[[1]]@DataMatrix
      }else{
        output <- matrix(NA, nrow=0, ncol=0)
      }
    }
  }else if(identical(typematch, "mRNA_Array")){
    if(is(object@mRNAArray, "list")){
      if(length(object@mRNAArray) > 1){
        output <- lapply(object@mRNAArray, function(tmp){getElement(tmp, "DataMatrix")})
        keeplist <- which.max(sapply(output, ncol))
        output <- output[[keeplist]]
        warning(paste("Taking the mRNA array platform with the greatest number of samples:", keeplist))
      }else if(length(object@mRNAArray)==1){
        output <- object@mRNAArray[[1]]@DataMatrix
      }else{
        output <- matrix(NA, nrow=0, ncol=0)
      }
    }
  }else if(identical(typematch, "GISTIC_A")){
    output <- getElement(object@GISTIC, "AllByGene")
  }else if(identical(typematch, "GISTIC_T")){
    output <- getElement(object@GISTIC, "ThresholedByGene")
  }else{
    stop(paste("Type", typematch, "not yet supported."))
  }
  if(dim(output)[1] == 0 | dim(output)[2] == 0){
    message("There is no data for that data type!")
  } else {
    
    if (phenoData){
    bcID <- function(barcodes, center=FALSE){
      if(center){
        bcc <- sapply(strsplit(barcodes, "-"), "[", c(4:7))
        bcc <- tolower(t(bcc))
      } else {
        bcc <- lapply(strsplit(barcodes, "-"), "[", c(1:3))
        bcc <- tolower(sapply(bcc, paste, collapse="-"))
        bcc <- as.vector(bcc)
      }
      return(bcc)
    }
    
    dm <- apply(output[4:ncol(output)], 2, as.numeric, as.matrix)
    rownames(dm) <- rownames(output)
    dups <- bcID(colnames(dm))[duplicated(bcID(colnames(dm)))]
    
    d <- c()
    for (cc in 1:length(dups)){
      d <- cbind(d, apply(dm[,colnames(dm) %in% dups[cc]], 1, mean))
    }
    colnames(d) <- dups
    dm <- cbind(dm[,!(bcID(colnames(dm)) %in% dups)], d)
    
    centerbc <- bcID(colnames(output)[-c(1:3)], center=TRUE)
    colnames(centerbc) <- c("sample", "portion", "plate", "center")
    rownames(centerbc) <- bcID(colnames(output))[-c(1:3)]
    
    pd <- getElement(object, "Clinical")
    if(length(pd)==0){
      stop("No clinical data available!")
    }
    
    npd <- pd[na.omit(match(colnames(dm), rownames(pd))),]
    npd <- cbind(npd, centerbc[match(rownames(npd), rownames(centerbc)),])
    
    ndm <- dm[,na.omit(match(rownames(pd), colnames(dm)))]
    
    if(identical(all.equal(rownames(npd), colnames(ndm)), TRUE)){
      eset <- ExpressionSet(ndm, AnnotatedDataFrame(npd))
    }else{
      stop("Couldn't match up rownames of phenodata to colnames of numeric data")
    }
    
    
      ovextraclin <- fread("extraclinfilehere")
      ovextraclin[grep("patient_barcode", ovextraclin[, 1]),][-1]
      colnames(ovextraclin)[-1] <- ovextraclin[grep("patient_barcode", ovextraclin[, 1]),][-1]
      rownames(ovextraclin) <- ovextraclin[, 1]
      
      jj <- ovextraclin[,match(rownames(npd), colnames(ovextraclin))]
      featureData(eset) <- AnnotatedDataFrame(jj)
      return(eset)
    } else {
      return(output)
    }
  }
}