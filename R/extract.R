.dlClinx <- function(object){
  adt <- getElement(object, "runDate")
  dset <- getElement(object, "Dataset")
  cl_url <- "http://gdac.broadinstitute.org/runs/stddata__"
  cl_url <- paste0(cl_url,substr(adt,1,4),"_",substr(adt,5,6),"_",substr(adt,7,8),"/data/")
  cl_url <- paste0(cl_url,dset,"/",adt,"/")
  cl_url <- paste0(cl_url, "gdac.broadinstitute.org_", dset, ".Merge_Clinical.Level_1.", adt, "00.0.0.tar.gz")
  
  download.file(url=cl_url, destfile=paste0(dset, "-ExClinical.tar.gz"), method="auto", quiet=TRUE, mode="w")
  fileList <- untar(paste0(dset, "-ExClinical.tar.gz"), list=TRUE)
  fileList <- grep(".clin.merged.txt", fileList, fixed = TRUE, value=TRUE)
  untar(paste0(dset,"-ExClinical.tar.gz"),files=fileList)
  filename <- paste0(adt,"-",dset,"-ExClinical.txt")
  file.rename(from=fileList,to=filename)
  file.remove(paste0(dset,"-ExClinical.tar.gz"))
  unlink(strsplit(fileList[1],"/")[[1]][1], recursive = TRUE)
  extracl <- data.table::fread(filename, data.table=FALSE)
  file.remove(filename)
  colnames(extracl)[-1] <- extracl[grep("patient_barcode", extracl[, 1]),][-1]
  rownames(extracl) <- extracl[, 1]      
  extracl <- extracl[,-1]
  extracl <- t(extracl)
  extracl <- extracl[,!grepl("patient_barcode", colnames(extracl))]
  return(extracl)
}

.fileSelect <- function() {
  g <- readline(
    paste0("The selected data type has more than one", 
           "file available.\nPlease select the desired file.",
           "\n(Enter 0 for the first file with the most number of samples)\n_"))
  g <- suppressWarnings(as.integer(g))
  if(is.na(g)){
    stop("Your selection must be an integer!")
  } else {
    return(g) 
  }}
#' Extract data from \code{FirehoseData} object into \code{ExpressionSet} or \code{GRangesList} object
#' 
#' This extracts data from a \code{FirehoseData} object and converts it to a structured S4 
#' object for analysis. An option is available to retreive clinical data. The function 
#' returns an ExpressionSet, GRangesList class object (see \code{\link{ExpressionSet}} or \code{\link{GRangesList}})
#' or a data frame. 
#' 
#' @param object A \code{FirehoseData} object from which to extract data. 
#' @param type The type of data to extract from the "FirehoseData" object. To request the clinical data frame only, set to NULL or "none".
#' @param clinical Logical (default TRUE) includes additional clinic data, if available.
#' @return Either an \code{\link{ExpressionSet}} object, \code{\link{GRangesList}}) object, or data frame (when no clinical data requested for the selected data type or when clinical data is exclusively requested). 
#' Choices include: "RNAseq_Gene", "Clinic", "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA", "GISTIC_A", "GISTIC_T".
#' The "GISTIC_A" type of dataset represents GISTIC data by all genes. "GISTIC_T"" represents data thresholded by genes.
#' 
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#' 
#' @examples 
#' 
#' \dontrun{
#' b2 <- extract(a2, "Methylation", clinical=TRUE)
#' }
#' 
#' @export
extract <- function(object, type, clinical = TRUE){
  if (!is.null(type)) {
    if (is.character(type)) {
      type <- tolower(gsub("_", "", type)) 
    } else {
      stop("Data type must be a character string")
    }
  } else {
    if (!clinical){
      stop("Nothing to extract. Please check the arguments.")
    }
  }
  if(clinical){
    pd <- getElement(object, "Clinical")
    if(any(grepl("\\.", rownames(pd)))){
      rownames(pd) <- gsub("\\.", "-", rownames(pd))
    }
    if(length(pd)==0){
      stop("No clinical data available!")
    }
    clinextra <- .dlClinx(object)
    if(any(grepl("\\.", rownames(clinextra)))){
      rownames(clinextra) <- gsub("\\.", "-", rownames(clinextra))
    }
    if(length(clinextra)==0){
      message("No additional clinical information found!")
    }
    pd <- merge(pd, clinextra, "row.names")
    rownames(pd) <- pd[,"Row.names"]
    pd <- pd[,-c(1,2)]
    if(is.null(type) || type == "none") {	
      return(pd) 
    }
  } 
  if (grepl("s$", type)) {
    gsub("s$", "", type)
  }
  choices <- tolower(
    gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
                    "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
                    "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
                    "miRNA_Array", "RPPA")))
  if(type %in% choices){
    slotreq <- grep(paste0("^", type) , slotNames(object), 
                    ignore.case=TRUE, perl=TRUE, value=TRUE)
    if(is(getElement(object, slotreq), "list")){
      elemlength <- length(getElement(object, slotreq))
      if(elemlength > 1){
        if(interactive()){
          sourceName <- sapply(getElement(object, slotreq), function(FHarray) { getElement(FHarray, "Filename") } )
          dimensions <- sapply(lapply(getElement(object, slotreq), function(tmp){getElement(tmp, "DataMatrix")}), dim)
          cat(paste0("[", seq(length(sourceName)), "] ", sourceName, paste0("\n\tNumber of rows: ", dimensions[1,], "\tNumber of columns: ", dimensions[2,]) ), fill = TRUE, sep = "\n")
          fileNo <- .fileSelect()
          if(fileNo == 0) { fileNo <- which.max(sapply(lapply(getElement(object, slotreq), function(tmp){getElement(tmp, "DataMatrix")}), ncol)) }
          message("Selecting file: [", fileNo, "] ", sourceName[fileNo])
          dm <- getElement(object, slotreq)[[fileNo]]@DataMatrix
        } else {
          dm <- lapply(getElement(object, slotreq), function(tmp){getElement(tmp, "DataMatrix")})
          keeplist <- which.max(sapply(dm, ncol))
          dm <- dm[[keeplist]]
          warning(paste("Taking the array platform with the greatest number of samples:", keeplist))
        }
      } else if(elemlength == 1){
        dm <- getElement(object, slotreq)[[1]]@DataMatrix
      } else if(elemlength == 0){
        dm <- matrix(NA, nrow=0, ncol=0)
      }
    } else {
      dm <- getElement(object, slotreq)
    }
  } else if(type %in% c("gistica", "gistict")) {
    if(type=="gistica"){
      slotreq <- "AllByGene"
    } else {
      slotreq <- "ThresholdedByGene"
    }
    dm <- getElement(object@GISTIC, slotreq)
  } else {
    stop(paste("Data type not yet supported or could not be matched."))
  }
  if(dim(dm)[1] == 0 | dim(dm)[2] == 0){
    stop("There is no data for that data type!")
  } else {
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH", "Mutations")
    if(slotreq %in% c("Methylation", "AllByGene", "ThresholdedByGene")){
      annote <- dm[, seq(grep("TCGA", names(dm))[1]-1) ]
      rNames <- rownames(dm)
      dm <- apply(dm[grep("TCGA", names(dm))[1]:ncol(dm)], 2, as.numeric, as.matrix)
      rownames(dm) <- rNames
      if(any(grepl("\\.", colnames(dm)))){ colnames(dm) <- gsub("\\.", "-", colnames(dm)) }
      righttab <- bcRight(colnames(dm))
      dups <- bcIDR(colnames(dm), sample=T, collapse = T)[duplicated(bcIDR(colnames(dm), sample = T, collapse = T))]
    } else if(slotreq %in% rangeslots) {
      if(slotreq=="Mutations"){
        dm <- split(dm, dm$Tumor_Sample_Barcode)
      } else {
        dm <- split(dm, dm$Sample)
      } 
    ## Checking for duplicates in data
      if(any(grepl("\\.", names(dm)))){ names(dm) <- gsub("\\.", "-", names(dm)) }
      dups <- bcIDR(names(dm), sample=T, collapse=T)[duplicated(bcIDR(names(dm), sample=T, collapse = T))]
    } else {
      dups <- bcIDR(colnames(dm), sample=T, collapse=T)[duplicated(bcIDR(colnames(dm), sample=T, collapse = T))]
    }
    # check for presence of technical replicates
    if(length(dups) != 0){ 
      if(!slotreq %in% rangeslots){
        repeated <- dm[, bcIDR(colnames(dm), sample=T, collapse=T) %in% dups]
        duplic <- bcIDR(colnames(repeated), sample=T, collapse=T)[duplicated(bcIDR(colnames(repeated), sample=T, collapse=T))]
        d <- c()
        for (cc in seq(duplic)){
          d <- cbind(d, apply(repated[,bcIDR(colnames(repeated)) %in% bcIDR(duplic[cc])], 1, mean))
        }
        colnames(d) <- bcIDR(duplic)
        dm <- cbind(dm[, !(bcIDR(colnames(dm)) %in% bcIDR(dups))], d)
        colnames(dm) <- bcIDR(colnames(dm))
      } else {
        repeated <- names(dm)[bcIDR(names(dm), sample=T, collapse=T) %in% dups]
        # dropping technical replicates 
        dm <- dm[!duplicated(bcIDR(names(dm), sample=T, collapse=T))]
        duplic <- bcIDR(repeated, sample=T, collapse=T)[duplicated(bcIDR(repeated, sample=T, collapse=T))]
      }
    }
    if(!slotreq %in% rangeslots){
      righttab <- bcRight(colnames(dm))
    } else {
      righttab <- bcRight(names(dm))
    }
    cleanDupCols <- function(object){
      if(is(object, "list")){
        cleanObj <- lapply(object, function(UBC) {
          UBC[!duplicated(lapply(UBC, c))]
        }) 
        return(cleanObj)
      } else if(is.data.frame(object)){
        object <- object[,!duplicated(lapply(object, c))]
        return(object)
      }
    }
    if(exists("pd")){
    pd <- cleanDupCols(pd)
    if(!slotreq %in% rangeslots){
      clindup <- matrix(NA, nrow=ncol(dm))
      rownames(clindup) <- bcIDR(colnames(dm), sample=TRUE, collapse=TRUE)
      clindup <- cbind(clindup, pd[match(bcIDR(rownames(clindup)), bcIDR(rownames(pd))),])
      rowsNotEmpty <- apply(clindup, 1, function(x) !all(is.na(x)))
      clindup <- clindup[rowsNotEmpty,]
      clindup <- clindup[, -(which(colnames(clindup) == "clindup"))]
      righttable <- righttab[match(rownames(clindup), righttab[, "patientids"]),]
      righttable <- righttable[, -(which(colnames(righttable) == "patientids"))]
      clindup <- cbind(clindup, righttable)
      colnames(dm) <- bcIDR(colnames(dm), sample=TRUE, collapse=TRUE)
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
    } else {
      clindup <- matrix(NA, nrow=length(dm))
      rnames <- bcIDR(names(dm), sample=TRUE, collapse=TRUE)
      clindup <- pd[match(bcIDR(rnames), bcIDR(rownames(pd))),]
      clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
      clindup <- cbind(clindup, righttab[match(rnames, rownames(righttab)),])
      clindup <- S4Vectors::DataFrame(patientids = rownames(clindup), clindup, row.names = NULL)
      names(dm) <- bcIDR(names(dm), sample=TRUE, collapse=TRUE)
      matchLogic <- bcIDR(names(dm)) %in% clindup[, "patientids"]
      if (all(!matchLogic)) {
        stop("no names could be matched")
      }
      dm <- dm[matchLogic]
      dm <- lapply(dm, FUN = function(ubc) { names(ubc) <- tolower(names(ubc)) 
      ubc } )
      if(slotreq=="Mutations"){
        metdat <- lapply(dm, FUN = function(lt) {
          lt <- lt[, !names(lt) %in% c("chromosome","start_position",
                                       "end_position", "strand", 
                                       "center","start","end","gene",
                                       "dataset","type","chr","ref_allele",
                                       "tum_allele1","tum_allele2")]
          return(lt)
        })
        mygrl <- GRangesList(lapply(dm, FUN = function(gr){
          GRanges(seqnames = paste0("chr", Rle(gr$chromosome)), 
                  ranges = IRanges(as.numeric(gr$start_position),
                                   as.numeric(gr$end_position),
                  names = gr$hugo_symbol),
                  strand = gr$strand)}
        ))
        ncbi_build <- as.character(Reduce(intersect, 
                                          lapply(metdat, function(x)
                                          {
                                            x[, "ncbi_build"]
                                          })))
        if (length(ncbi_build) == 1L) {
          genome(mygrl) <- ncbi_build
        } else {
          message("NCBI build was not consistent")
        }
        metdat <- cleanDupCols(metdat)
        mygrl@metadata <- metdat
      } else {
        mygrl <- GRangesList(lapply(dm, FUN = function(gr){
          GRanges(seqnames = paste0("chr", Rle(gr$chromosome)), 
                  ranges = IRanges(gr$start, gr$end),
                  Num_Probes = gr$num_probes, 
                  Segment_Mean = gr$segment_mean)}
        ))
      }
      mcols(mygrl) <- clindup
      if(exists("sourceName")) {
        mygrl@metadata <- list("fileName" = sourceName[fileNo])
      }
      return(mygrl)
    }
    } else {
    return(dm)
    }
  }
}