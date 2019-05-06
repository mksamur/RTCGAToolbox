.getListData <- function(object, platform) {
    if (is.null(platform)) { stop("Please set platform") }
    switch(class(platform),
        "numeric" = {
            if (platform > length(object)) {
                message("Accessible platforms:")
                for (i in seq_along(object)) {
                    message(paste0("#",i," :",object[[i]]@Filename))
                }
                stop("Invalid list member")
            }
            invisible(object[[platform]]@DataMatrix)
        }, {
            message("Accessible platforms:")
            for (i in seq_along(object)) {
                message(paste0("#",i," :",object[[i]]@Filename))
            }
            stop("Set a valid platfrom")
        }
    )
}

#' An S4 class to store data from CGA platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data frame that stores the CGH data.
#' @exportClass FirehoseCGHArray
setClass("FirehoseCGHArray", representation(Filename = "character",
    DataMatrix = "data.frame"))
setMethod("show", "FirehoseCGHArray",function(object) {
    message(paste0("Platform:", object@Filename))
    if (dim(object@DataMatrix)[1] > 0 ) {
        message("FirehoseCGHArray object, dim: ", paste(dim(object@DataMatrix),
        collapse = "\t"))
    }
})

#' An S4 class to store data from methylation platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data frame that stores the methylation data.
#' @exportClass FirehoseMethylationArray
setClass("FirehoseMethylationArray", representation(Filename = "character", DataMatrix = "data.frame"))
setMethod("show", "FirehoseMethylationArray",function(object){
  message(paste0("Platform:", object@Filename))
  if(dim(object@DataMatrix)[1] > 0 ){message("FirehoseMethylationArray object, dim: ", paste(dim(object@DataMatrix),collapse = "\t"))}
})


#' An S4 class to store data from array (mRNA, miRNA etc.) platforms
#'
#' @slot Filename Platform name
#' @slot DataMatrix A data matrix that stores the expression data.
#' @exportClass FirehosemRNAArray
setClass("FirehosemRNAArray", representation(Filename = "character", DataMatrix = "matrix"))
setMethod("show", "FirehosemRNAArray",function(object){
  message(object@Filename)
  if(dim(object@DataMatrix)[1] > 0 ){message("FirehoseCGHArray object, dim: ", paste(dim(object@DataMatrix),collapse = "\t"))}
})

#' An S4 class to store processed copy number data. (Data processed by using GISTIC2 algorithm)
#'
#' @slot Dataset Cohort name
#' @slot AllByGene A data frame that stores continuous copy number
#' @slot ThresholdedByGene A data frame for discrete copy number data
#' @slot Peaks A data frame storing GISTIC peak data.
#' See \link{getGISTICPeaks}.
#' @exportClass FirehoseGISTIC
setClass("FirehoseGISTIC", representation(Dataset = "character",
    AllByGene = "data.frame", ThresholdedByGene="data.frame",
    Peaks = "data.frame"))

setMethod("show", "FirehoseGISTIC", function(object){
    if (.hasOldGISTIC(object)) {
        warning("'FirehoseGISTIC' object is outdated, please run 'updateObject()'")
    }
    message(paste0("Dataset:", object@Dataset))
    if (dim(object@AllByGene)[1L] > 0L) {
        message("FirehoseGISTIC object, dim: ", paste(dim(object@AllByGene),
            collapse = "\t"))
    }
})

#' @importFrom S4Vectors isEmpty
#' @importFrom methods slotNames
#' @param x A FirehoseGISTIC class object
#' @exportMethod isEmpty
#' @describeIn FirehoseGISTIC check whether the FirehoseGISTIC object has
#' data in it or not
setMethod("isEmpty", "FirehoseGISTIC", function(x) {
    allSlots <- slotNames(x)
    all(vapply(allSlots, function(g) {
        obj <- getElement(x, g)
        if (is.data.frame(obj) || is.character(obj))
            !length(obj)
        else
            isEmpty(obj)
    }, logical(1L)))
})

#' An S4 class to store main data object from clinent function.
#'
#' @slot Dataset A cohort name
#' @slot runDate Standard data run date from \code{\link{getFirehoseRunningDates}}
#' @slot gistic2Date Analyze running date from \code{\link{getFirehoseAnalyzeDates}}
#' @slot clinical clinical data frame
#' @slot RNASeqGene Gene level expression data matrix from RNAseq
#' @slot RNASeq2GeneNorm Gene level expression data matrix from RNAseq (RSEM)
#' @slot miRNASeqGene miRNA expression data from matrix smallRNAseq
#' @slot CNASNP A data frame to store somatic copy number alterations from SNP array platform
#' @slot CNVSNP A data frame to store germline copy number variants from SNP array platform
#' @slot CNASeq A data frame to store somatic copy number alterations from sequencing platform
#' @slot CNACGH A list that stores \code{FirehoseCGHArray} object for somatic copy number alterations from CGH platform
#' @slot Methylation A list that stores \code{FirehoseMethylationArray} object for methylation data
#' @slot mRNAArray A list that stores \code{FirehosemRNAArray} object for gene expression data from microarray
#' @slot miRNAArray A list that stores \code{FirehosemRNAArray} object for miRNA expression data from microarray
#' @slot RPPAArray A list that stores \code{FirehosemRNAArray} object for RPPA data
#' @slot Mutation A data frame for mutation infromation from sequencing data
#' @slot GISTIC A \code{FirehoseGISTIC} object to store processed copy number data
#' @slot BarcodeUUID A data frame that stores the Barcodes, UUIDs and Short sample identifiers
#' @exportClass FirehoseData
setClass("FirehoseData", representation(Dataset = "character",
    runDate = "character", gistic2Date = "character", clinical = "data.frame",
    RNASeqGene = "matrix", RNASeq2GeneNorm="list", miRNASeqGene="matrix",
    CNASNP="data.frame", CNVSNP="data.frame", CNASeq="data.frame", CNACGH="list",
    Methylation="list", mRNAArray="list", miRNAArray="list", RPPAArray="list",
    Mutation="data.frame", GISTIC="FirehoseGISTIC", BarcodeUUID="data.frame"))

#' @describeIn FirehoseData show method
#' @param object A FirehoseData object
setMethod("show", "FirehoseData",function(object) {
    if (.hasOldAPI(object) || .hasOldGISTIC(selectType(object, "GISTIC"))) {
        object <- updateObject(object)
    warning("'FirehoseData' object is outdated, please run 'updateObject()'")
    }
    cat(paste0(object@Dataset," FirehoseData object"))
    cat(paste0("Standard run date: ", object@runDate), "\n")
    cat(paste0("Analysis running date: ", object@gistic2Date), "\n")
    cat("Available data types:", "\n")
    if (dim(object@clinical)[1] > 0 & dim(object@clinical)[2] > 0) {
        cat("  clinical: A data frame of phenotype data, dim: ",
            paste(dim(object@clinical), collapse = " x "), "\n")}
    if (dim(object@RNASeqGene)[1] > 0 & dim(object@RNASeqGene)[2] > 0) {
        cat("  RNASeqGene: A matrix of count or normalized data, dim: ",
            paste(dim(object@RNASeqGene),collapse = " x "), "\n")}
    if (length(object@RNASeq2GeneNorm)) {
        cat("  RNASeq2GeneNorm: A list of FirehosemRNAArray object(s), length: ",
            length(object@RNASeq2GeneNorm), "\n")}
    if (dim(object@miRNASeqGene)[1] > 0 & dim(object@miRNASeqGene)[2] > 0) {
        cat("  miRNASeqGene: A matrix, dim: ",
            paste(dim(object@miRNASeqGene), collapse = " x "), "\n")}
    if (dim(object@CNASNP)[1] & dim(object@CNASNP)[2]) {
        cat("  CNASNP: A data.frame, dim: ", paste(dim(object@CNASNP),
            collapse = " x "), "\n")}
    if (dim(object@CNVSNP)[1] & dim(object@CNVSNP)[2]) {
        cat("  CNVSNP: A data.frame, dim: ", paste(dim(object@CNVSNP),collapse = " x "), "\n")}
    if (dim(object@CNASeq)[1] & dim(object@CNASeq)[2]) {
        cat("  CNASeq: A data.frame, dim: ", paste(dim(object@CNASeq),collapse = " x "), "\n")}
    if (length(object@CNACGH)) {
        cat("  CNACGH: A list of FirehoseCGHArray object(s), length: ",
                length(object@CNACGH), "\n")}
    if (length(object@Methylation)) {
        cat("  Methylation: A list of FirehoseMethylationArray object(s), length: ",
                length(object@Methylation), "\n")}
    if (length(object@mRNAArray)) {
        cat("  mRNAArray: A list of FirehosemRNAArray object(s), length: ",
                length(object@mRNAArray), "\n")}
    if (length(object@miRNAArray)) {
        cat("  miRNAArray: A list of FirehosemRNAArray object(s), length: ",
                length(object@miRNAArray), "\n")}
    if (length(object@RPPAArray)) {
        cat("  RPPAArray: A list of FirehosemRNAArray object(s), length: ",
                length(object@RPPAArray), "\n")}
    if (length(object@GISTIC@Dataset)) {
        cat("  GISTIC: A FirehoseGISTIC for copy number data", "\n")}
    if (dim(object@Mutation)[2] & dim(object@Mutation)[2]) {
        cat("  Mutation: A data.frame, dim: ", paste(dim(object@Mutation),
        collapse = " x "), "\n")}
    cat("To export data, use the 'getData' function.\n")
})

#' @title Extract data from FirehoseData object
#'
#' @description A go-to function for getting top level information from a
#' \code{\linkS4class{FirehoseData}} object. Available datatypes for a
#' particular object can be seen by entering the object name in the
#' console ('show' method).
#'
#' @param object A \code{\linkS4class{FirehoseData}} object
#' @param type A data type to be extracted
#' @param platform An index for data types that may come from multiple
#' platforms (such as mRNAArray), for GISTIC data, one of the options:
#' 'AllByGene' or 'ThresholdedByGene'
#'
#' @examples
#' data(RTCGASample)
#' getData(RTCGASample, "clinical")
#' getData(RTCGASample, "RNASeqGene")
#'
#' @return Returns matrix or data.frame depending on data type
setGeneric("getData", function(object, type, platform) {
    standardGeneric("getData")
})

#' @describeIn FirehoseData Get a matrix or data.frame from \code{FirehoseData}
#' @param type A data type to be extracted
#' @param platform An index for data types that may come from multiple
#' platforms (such as mRNAArray), for GISTIC data, one of the options:
#' 'AllByGene' or 'ThresholdedByGene'
#' @exportMethod getData
setMethod("getData", "FirehoseData",
    function(object, type, platform) {
        withPlat <- c("CNACGH", "mRNAArray", "Methylation", "miRNAArray",
            "RPPAArray")
        if (missing(type))
            stop("Enter a data type")
        if (type %in% withPlat) {
            res <- .getListData(getElement(object, type), platform)
        } else if (identical(type, "GISTIC")) {
            if (!platform %in% c("ThresholdedByGene", "AllByGene", "Peaks") ||
                !S4Vectors::isSingleString(platform))
        stop("GISTIC platforms available:\n",
             "\t'AllByGene', 'ThresholdedByGene', & 'Peaks'")
            res <-  getElement(selectType(object, "GISTIC"), platform)
        } else {
            res <- selectType(object, type)
        }
        if (!length(res))
            stop("No data available for that type")
        res
})

#' An S4 class to store differential gene expression results
#'
#' @slot Dataset Dataset name
#' @slot Toptable Results data frame
#' @exportClass DGEResult
setClass("DGEResult", representation(Dataset = "character", Toptable = "data.frame"))
setMethod("show", "DGEResult",function(object){
  message(paste0("Dataset:", object@Dataset))
  if(dim(object@Toptable)[1] > 0 ){message("DGEResult object, dim: ", paste(dim(object@Toptable),collapse = "\t"))}
})

#' Export toptable or correlation data frame
#' @param object A \code{\linkS4class{DGEResult}} or \code{\linkS4class{CorResult}} object
#' @return Returns toptable or correlation data frame
#' @examples
#' data(RTCGASample)
#' dgeRes = getDiffExpressedGenes(RTCGASample)
#' dgeRes
#' showResults(dgeRes[[1]])
setGeneric("showResults",
           function(object) standardGeneric("showResults")
)

#' Export toptable or correlation data frame
#' @param object A \code{\linkS4class{DGEResult}} or \code{\linkS4class{CorResult}} object
#' @rdname showResults-DGEResult
#' @aliases showResults,DGEResult,DGEResult-method
#' @return Returns toptable for DGE results
#' @export
#' @examples
#' data(RTCGASample)
#' dgeRes = getDiffExpressedGenes(RTCGASample)
#' dgeRes
#' showResults(dgeRes[[1]])
setMethod("showResults", "DGEResult",function(object){
  message(paste0("Dataset: ",object@Dataset))
  print(head(object@Toptable))
  invisible(object@Toptable)
})

#' An S4 class to store correlations between gene expression level and copy number data
#'
#' @slot Dataset A cohort name
#' @slot Correlations Results data frame
#' @exportClass CorResult
setClass("CorResult", representation(Dataset = "character", Correlations = "data.frame"))
setMethod("show", "CorResult",function(object){
  message(paste0("Dataset:", object@Dataset))
  if(dim(object@Correlations)[1] > 0 ){message("CorResult object, dim: ", paste(dim(object@Correlations),collapse = "\t"))}
})

#' Export toptable or correlation data frame
#' @param object A \code{\linkS4class{DGEResult}} or \code{\linkS4class{CorResult}} object
#' @rdname showResults-CorResult
#' @aliases showResults,CorResult,CorResult-method
#' @return Returns correlation results data frame
#' @examples
#' data(RTCGASample)
#' corRes = getCNGECorrelation(RTCGASample,adj.pval = 1,raw.pval = 1)
#' corRes
#' showResults(corRes[[1]])
setMethod("showResults", "CorResult",function(object){
  message(paste0("Dataset: ",object@Dataset))
  print(head(object@Correlations))
  invisible(object@Correlations)
})


.hasOldAPI <- function(object) {
    isTRUE(methods::.hasSlot(object, "RNAseq")) ||
    isTRUE(methods::.hasSlot(object, "Mutations")) ||
    isTRUE(methods::.hasSlot(object, "Clinical"))
}

.hasOldGISTIC <- function(object) {
    isTRUE(methods::.hasSlot(object, "ThresholedByGene"))
}

#' @describeIn FirehoseData Update an old RTCGAToolbox FirehoseData object to
#' the most recent API
#' @param verbose logical (default FALSE) whether to print extra messages
#' @param ... additional arguments for updateObject
#' @exportMethod updateObject
setMethod("updateObject", "FirehoseData",
    function(object, ..., verbose = FALSE) {
    if (verbose)
        message("updateObject(object = 'FirehoseData')")
    oldAPI <- .hasOldAPI(object)
    oldGISTIC <- .hasOldGISTIC(selectType(object, "GISTIC"))
    if (oldAPI) {
    object <- new(class(object), Dataset = object@Dataset,
        runDate = NA_character_, gistic2Date = NA_character_,
        clinical = if (.hasSlot(object, "Clinical")) { object@Clinical }
        else { object@clinical },
        RNASeqGene = object@RNASeqGene,
        RNASeq2GeneNorm = object@RNASeq2GeneNorm,
        miRNASeqGene = object@miRNASeqGene, CNASNP = object@CNASNP,
        CNVSNP = object@CNVSNP,
        CNASeq = if (.hasSlot(object, "CNAseq")) { object@CNAseq }
        else { object@CNASeq },
        CNACGH = object@CNACGH, Methylation = object@Methylation,
        mRNAArray = object@mRNAArray, miRNAArray = object@miRNAArray,
        RPPAArray = object@RPPAArray,
        Mutation = if (.hasSlot(object, "Mutations")) { object@Mutations }
        else { object@Mutation },
        GISTIC = object@GISTIC, BarcodeUUID = object@BarcodeUUID)
    }
    if (oldGISTIC) {
       object@GISTIC <- updateObject(selectType(object, "GISTIC"))
    }
    return(object)
})

#' @describeIn FirehoseGISTIC Update an old FirehoseGISTIC object to
#' the most recent API
#' @param object A \code{FirehoseGISTIC} object
#' @param verbose logical (default FALSE) whether to print extra messages
#' @param ... additional arguments for updateObject
#' @exportMethod updateObject
setMethod("updateObject", "FirehoseGISTIC",
    function(object, ..., verbose = FALSE) {
    if (verbose)
        message("updateObject(object = 'FirehoseGISTIC')")
    oldGISTIC <- .hasOldGISTIC(object)
    if (oldGISTIC) {
        object <- new("FirehoseGISTIC", Dataset = object@Dataset,
            AllByGene = object@AllByGene,
            ThresholdedByGene = object@ThresholedByGene)
    }
    return(object)
})
