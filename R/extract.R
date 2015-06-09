#' Extract data from \code{FirehoseData} object into S4 \code{ExpressionSet}
#' 
#' This function serves to extract and reorganize data into a structured S4 
#' object for analysis. An option is available to retreive clinical data 
#' and return an ExpressionSet class object (see \code{\link{ExpressionSet}})
#' otherwise, a data frame is returned. 
#' 
#' @param object A \code{FirehoseData} object from which to extract data. 
#' @param type The type of data to extract from the "FirehoseData" object.
#' @param clinical Logical (default TRUE) includes additional clinic data, if available.
#' @return An \code{\link{ExpressionSet}} object for the selected data type or data frame if no clinical data is requested. 
#' Choices include: "RNAseq_Gene", "Clinic", "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA", "GISTIC_A", "GISTIC_T".
#' The "GISTIC_A" type of dataset represents GISTIC data by all genes. "GISTIC_T"" represents data thresholded by genes.
#' 
#' @examples 
#' 
#' \dontrun{
#' b2 <- extract(a2, "Methylation", clinical=TRUE)
#' }
#' 
extract <- function(object, type, clinical = TRUE){
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
      slotreq <- "AllByGene"
    } else {
      slotreq <- "TresholedByGene"
    }
    dm <- getElement(object@GISTIC, slotreq)
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
      } else if(sample & !collapse){
        samp <- tolower(sapply(strsplit(barcodes, "-"), "[", 4))
        samp <- as.numeric(substr(samp, 1, nchar(samp)-1))
        return(samp)
      }
    }
    return(bcc)
  }
  # require just character vector input
  rightbc <- function(mat){
    dmat <- dm
    if(slotreq %in% rangeslots && is(dm, "list")){
      FUN = names
    } else {
      FUN = colnames
    }
    sample_type <- samptab[,2][match(bcID(FUN(dmat), sample=TRUE), samptab[,1])]
    tb <- data.frame(sample_type, sample_code = substr(bcID(FUN(dmat), sample=TRUE), 1,2), vial = substr(bcID(FUN(dmat), sample=TRUE, center=TRUE)[,1], 3,3),
                     portion = substr(bcID(FUN(dmat), portion=TRUE), 1,2), analyte = substr(bcID(FUN(dmat), portion=TRUE), 3,3), 
                     plate = bcID(FUN(dmat), center=T)[,1], center = bcID(FUN(dmat), center=T)[,2],  stringsAsFactors=FALSE)
    tb[, "sample_code"] <- as.numeric(tb[, "sample_code"])
    rownames(tb) <- bcID(FUN(dmat), sample=TRUE, collapse=TRUE)
    return(tb)
  }
  
  if(dim(dm)[1] == 0 | dim(dm)[2] == 0){
    stop("There is no data for that data type!")
  } else {
    if(slotreq %in% c("Methylation", "AllByGene", "ThresholedByGene" )){
      annote <- dm[, seq(grep("TCGA", names(dm))[1]-1) ]
      rNames <- rownames(dm)
      dm <- apply(dm[grep("TCGA", names(dm))[1]:ncol(dm)], 2, as.numeric, as.matrix)
      rownames(dm) <- rNames
      if(any(grepl("\\.", colnames(dm)))){ colnames(dm) <- gsub("\\.", "-", colnames(dm)) }
      
      dups <- bcID(colnames(dm), sample=T, collapse = T)[duplicated(bcID(colnames(dm), sample = T, collapse = T))]
      
    }
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH")
    if(slotreq %in% rangeslots){
      if(any(grepl("\\.", dm[, "Sample"]))){ dm[, "Sample"] <- gsub("\\.", "-", dm[, "Sample"]) }
      dm <- split(dm, dm$Sample)
      sample_type <- samptab[,2][match(bcID(names(dm), sample=TRUE), samptab[,1])]
      
      dups <- names(dm)[duplicated(bcID(names(dm), sample=T,collapse=T))]
      
    }
    
    righttab <- rightbc(dm)
    # check if technical replicates and tumor/normals present
    if(length(dups) != 0){    
      if(!slotreq %in% rangeslots){
        repeated <- dm[, bcID(colnames(dm), sample=T, collapse=T) %in% dups]
        
        if(range(bcID(colnames(dm), sample=TRUE))[2] > 19){
          controls <- dm[, righttab$sample_code >= 20 & righttab$sample_code <= 29]
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
        dm <- cbind(dm[, !(bcID(colnames(dm)) %in% bcID(dups))], d)
        colnames(dm) <- bcID(colnames(dm))
      } else {
        replicates <- names(tumors)[gsub(".$", "", bcID(names(tumors), sample=T, collapse=T)) %in% gsub(".$", "", bcID(dups, sample=T, collapse=T))]
        # which(bcID(names(tumors), sample=T) %in% c(0:9))
        duplic <- bcID(replicates)[duplicated(bcID(replicates))]
#         d <- c()
#         for (cc in seq(duplic)){
#           d <- cbind(d, apply(replicates[bcID(rFun(replicates)) %in% bcID(duplic[cc])], 1, mean))
#         }
      }
    }
    
    adt <- getElement(object, "runDate")
    dset <- getElement(object, "Dataset")
    
    dlClinx <- function(object){
      
      cl_url <- "http://gdac.broadinstitute.org/runs/stddata__"
      cl_url <- paste0(cl_url,substr(adt,1,4),"_",substr(adt,5,6),"_",substr(adt,7,8),"/data/")
      cl_url <- paste0(cl_url,dset,"/",adt,"/")
      cl_url <- paste0(cl_url, "gdac.broadinstitute.org_", dset, ".Merge_Clinical.Level_1.", adt, "00.0.0.tar.gz")
      
      download.file(url=cl_url, destfile=paste0(dset, "-ExClinical.tar.gz"), method="auto", quiet=TRUE, mode="w")
      fileList <- untar(paste0(dset, "-ExClinical.tar.gz"), list=TRUE)
      fileList <- grep(".clin.merged.txt", fileList, fixed = TRUE, value=TRUE)
      untar(paste0(dset,"-ExClinical.tar.gz"),files=fileList)
      file.rename(from=fileList,to=paste0(adt,"-",dset,"-ExClinical.txt"))
      file.remove(paste0(dset,"-ExClinical.tar.gz"))
      unlink(strsplit(fileList[1],"/")[[1]][1], recursive = TRUE)
      extracl <- fread(paste0(adt,"-",dset,"-ExClinical.txt"), data.table=FALSE)
      colnames(extracl)[-1] <- extracl[grep("patient_barcode", extracl[, 1]),][-1]
      rownames(extracl) <- extracl[, 1]      
      extracl <- extracl[,-1]
      extracl <- t(extracl)
      extracl <- extracl[,!grepl("patient_barcode", colnames(extracl))]
      return(extracl)
    }
    if(!clinical){
      return(dm)
    } else if(clinical){
      pd <- getElement(object, "Clinical")
      if(any(grepl("\\.", rownames(pd)))){ rownames(pd) <- gsub("\\.", "-", rownames(pd)) }
      if(length(pd)==0){ stop("No clinical data available!") }
      clinextra <- dlClinx(object)
      if(any(grepl("\\.", rownames(clinextra)))){ rownames(clinextra) <- gsub("\\.", "-", rownames(clinextra)) }
      if(length(clinextra)==0){ message("No additional clinical information found!")}
      pd <- merge(pd, clinextra, "row.names")
      rownames(pd) <- pd[,"Row.names"]
      pd <- pd[,-c(1,2)]
      
      
      if(!slotreq %in% rangeslots){
        clindup <- matrix(NA, nrow=ncol(dm))
        rownames(clindup) <- bcID(colnames(dm), sample=T, collapse=T)
        clindup <- cbind(clindup, pd[match(bcID(rownames(clindup)), bcID(rownames(pd))),])
        clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
        clindup <- clindup[, -c(1:2)]
        clindup <- cbind(clindup, righttab[match(rownames(clindup), rownames(righttab)),])
        colnames(dm) <- bcID(colnames(dm), sample=T, collapse=T)
        
        ndm <- dm[,na.omit(match(rownames(clindup), colnames(dm)))]
        
        if(identical(all.equal(rownames(clindup), colnames(ndm)), TRUE)){
          eset <- ExpressionSet(ndm, AnnotatedDataFrame(clindup))
          if(exists("annote")){
            featureData(eset) <- AnnotatedDataFrame(annote)
          }
          return(eset)
        } else {
          stop("Couldn't match up rownames of clinical to colnames of numeric data")
          return(dm)
        }
      } else if(slotreq %in% rangeslots) {
        clindup <- matrix(NA, nrow=length(dm))
        rownames(clindup) <- bcID(names(dm), sample=T, collapse=T)
        clindup <- cbind(clindup, pd[match(bcID(rownames(clindup)), bcID(rownames(pd))),])
        clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
        clindup <- clindup[, -c(1:2)]
        ## rightbc should work for GRLs
        clindup <- cbind(clindup, righttab[match(rownames(clindup), rownames(righttab)),])
        names(dm) <- bcID(names(dm), sample=T, collapse=T)
        
        mygrl <- GRangesList(lapply(dm, FUN = function(gr){ GRanges(seqnames = paste0("chr", Rle(gr$Chromosome)), 
                                                                         ranges = IRanges(gr$Start, gr$End), Num_Probes = gr$Num_Probes, 
                                                                         Segment_Mean = gr$Segment_Mean)}) )
        mcols(mygrl) <- clindup
        return(mygrl)
        }
      }
    }
  } 
}

  