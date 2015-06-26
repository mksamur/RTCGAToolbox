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
#' @author Marcel Ramos
#' 
#' @examples 
#' 
#' \dontrun{
#' b2 <- extract(a2, "Methylation", clinical=TRUE)
#' }
#' 
#' @export
extract <- function(object, type, clinical = TRUE){
  if(!is.character(type)){
     stop("Data type must be a character string!")
  }
  type <- tolower(gsub("_", "", type))
  choices <- tolower(gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
               "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
               "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
               "miRNA_Array", "RPPA")))
  if(type %in% choices){
    slotreq <- grep(paste0("^", type) , slotNames(object), 
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
  } else if(type %in% c("gistica", "gistict")) {
    if(type=="gistica"){
      slotreq <- "AllByGene"
    } else {
      slotreq <- "TresholedByGene"
    }
    dm <- getElement(object@GISTIC, slotreq)
  } else {
    stop(paste("Data type not yet supported or could not be matched."))
  }
  
  if(dim(dm)[1] == 0 | dim(dm)[2] == 0){
    stop("There is no data for that data type!")
  } else {
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH", "Mutations")
    
    if(slotreq %in% c("Methylation", "AllByGene", "ThresholedByGene")){
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
      if(any(grepl("\\.", names(dm)))){ names(dm) <- gsub("\\.", "-", names(dm)) }
      dups <- bcIDR(names(dm), sample=T, collapse=T)[duplicated(bcIDR(names(dm), sample=T, collapse = T))]
      righttab <- bcRight(names(dm)) 
    } else {
      dups <- bcIDR(colnames(dm), sample=T, collapse=T)[duplicated(bcIDR(colnames(dm), sample=T, collapse = T))]
      righttab <- bcRight(colnames(dm))
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
        duplic <- bcIDR(repeated, sample=T, collapse=T)[duplicated(bcIDR(repeated, sample=T, collapse=T))]
        # what to do with technical replicates of range data types?
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
	pd <- cleanDupCols(pd)    
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
      } else {
        clindup <- matrix(NA, nrow=length(dm))
        rownames(clindup) <- bcIDR(names(dm), sample=T, collapse=T)
        clindup <- cbind(clindup, pd[match(bcIDR(rownames(clindup)), bcIDR(rownames(pd))),])
        clindup <- clindup[apply(clindup, 1, function(x) !all(is.na(x))),]
        clindup <- clindup[, -c(1:2)]
        clindup <- cbind(clindup, righttab[match(rownames(clindup), rownames(righttab)),])
        names(dm) <- bcIDR(names(dm), sample=T, collapse=T)
        dm <- dm[na.omit(match(rownames(clindup), names(dm)))]
	# dm <- cleanDupCols(dm)
 	dm <- lapply(dm, FUN = function(ubc) { names(ubc) <- tolower(names(ubc)) 
						ubc } )
        if(slotreq=="Mutations"){

	metdat <- lapply(dm, FUN = function(lt) { lt <- lt[, !names(lt) %in% c("chromosome","start_position", "end_position", "strand", 
		"center","start","end","gene","dataset","type","chr","ref_allele","tum_allele1","tum_allele2")]
							lt  } )
        mygrl <- GRangesList(lapply(dm, FUN = function(gr){ GRanges(seqnames = paste0("chr", Rle(gr$chromosome)), 
                                                                      ranges = IRanges(as.numeric(gr$start_position), as.numeric(gr$end_position)), strand = gr$strand 
                                                                      )}))
	metdat <- cleanDupCols(metdat)
	mygrl@metadata <- metdat	
        } else {
        mygrl <- GRangesList(lapply(dm, FUN = function(gr){ GRanges(seqnames = paste0("chr", Rle(gr$chromosome)), 
                                                                      ranges = IRanges(gr$start, gr$end), Num_Probes = gr$num_probes, 
                                                                      Segment_Mean = gr$segment_mean)} ))
        }
        mcols(mygrl) <- clindup
        return(mygrl)
      }
      }
    }
}
  
