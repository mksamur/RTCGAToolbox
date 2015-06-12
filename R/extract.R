#' Extract data from \code{FirehoseData} object into \code{ExpressionSet} or \code{GRagesList} object
#' 
#' This extracts data from a \code{FirehoseData} object and converts it to a structured S4 
#' object for analysis. An option is available to retreive clinical data. The function 
#' returns an ExpressionSet, GRangesList class object (see \code{\link{ExpressionSet}} or \code{\link{GRangesList}})
#' or a data frame. 
#' 
#' @param object A \code{FirehoseData} object from which to extract data. 
#' @param type The type of data to extract from the "FirehoseData" object.
#' @param clinical Logical (default TRUE) includes additional clinic data, if available.
#' @return Either an \code{\link{ExpressionSet}} object, \code{\link{GRangesList}}) object, or data frame (if no clinical data requested) for the selected data type. 
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
  
  if(dim(dm)[1] == 0 | dim(dm)[2] == 0){
    stop("There is no data for that data type!")
  } else {
    if(slotreq %in% c("Methylation", "AllByGene", "ThresholedByGene" )){
      annote <- dm[, seq(grep("TCGA", names(dm))[1]-1) ]
      rNames <- rownames(dm)
      dm <- apply(dm[grep("TCGA", names(dm))[1]:ncol(dm)], 2, as.numeric, as.matrix)
      rownames(dm) <- rNames
      if(any(grepl("\\.", colnames(dm)))){ colnames(dm) <- gsub("\\.", "-", colnames(dm)) }
      
      dups <- bcIDR(colnames(dm), sample=T, collapse = T)[duplicated(bcIDR(colnames(dm), sample = T, collapse = T))]
    }
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH")
    if(slotreq %in% rangeslots){
      if(any(grepl("\\.", dm[, "Sample"]))){ dm[, "Sample"] <- gsub("\\.", "-", dm[, "Sample"]) }
      dm <- split(dm, dm$Sample)
      sample_type <- samptab[,2][match(bcIDR(names(dm), sample=TRUE), samptab[,1])]
      
      dups <- names(dm)[duplicated(bcIDR(names(dm), sample=T,collapse=T))]
      righttab <- bcRight(names(dm)) 
    } else {
      righttab <- bcRight(colnames(dm))
    }

    # check if technical replicates and tumor/normals present
    if(length(dups) != 0){    
      if(!slotreq %in% rangeslots){
        repeated <- dm[, bcIDR(colnames(dm), sample=T, collapse=T) %in% dups]
        
        if(range(bcIDR(colnames(dm), sample=TRUE))[2] > 19){
          controls <- dm[, righttab$sample_code >= 20 & righttab$sample_code <= 29]
          replicates <- repeated[, !bcIDR(colnames(repeated)) %in% bcIDR(colnames(normals)) & 
                                   !bcIDR(colnames(repeated)) %in% bcIDR(colnames(controls))]
        } else {      
          replicates <- repeated[, !bcIDR(colnames(repeated)) %in% bcIDR(colnames(normals))]
        }
        duplic <- colnames(replicates)[duplicated(bcIDR(colnames(replicates)))]
        d <- c()
        for (cc in seq(duplic)){
          d <- cbind(d, apply(replicates[,bcIDR(colnames(replicates)) %in% bcIDR(duplic[cc])], 1, mean))
        }
        colnames(d) <- bcIDR(duplic)
        dm <- cbind(dm[, !(bcIDR(colnames(dm)) %in% bcIDR(dups))], d)
        colnames(dm) <- bcIDR(colnames(dm))
      } else {
        replicates <- names(dm)[gsub(".$", "", bcIDR(names(dm), sample=T, collapse=T)) %in% gsub(".$", "", bcIDR(dups, sample=T, collapse=T))]
        duplic <- bcIDR(replicates)[duplicated(bcIDR(replicates))]
#         d <- c()
#         for (cc in seq(duplic)){
#           d <- cbind(d, apply(replicates[bcIDR(rFun(replicates)) %in% bcIDR(duplic[cc])], 1, mean))
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
        rownames(clindup) <- bcIDR(colnames(dm), sample=T, collapse=T)
        clindup <- cbind(clindup, pd[match(bcIDR(rownames(clindup)), bcIDR(rownames(pd))),])
        clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
        clindup <- clindup[, -c(1:2)]
        clindup <- cbind(clindup, righttab[match(rownames(clindup), rownames(righttab)),])
        colnames(dm) <- bcIDR(colnames(dm), sample=T, collapse=T)
        
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
        rownames(clindup) <- bcIDR(names(dm), sample=T, collapse=T)
        clindup <- cbind(clindup, pd[match(bcIDR(rownames(clindup)), bcIDR(rownames(pd))),])
        clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
        clindup <- clindup[, -c(1:2)]
        clindup <- cbind(clindup, righttab[match(rownames(clindup), rownames(righttab)),])
        names(dm) <- bcIDR(names(dm), sample=T, collapse=T)
        dm <- dm[na.omit(match(rownames(clindup), names(dm)))]
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

  