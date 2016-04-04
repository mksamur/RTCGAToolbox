setClinical <- function (expData, phenodat) {
  if (is(expData, "list")) {
    names(expData) <- gsub("\\.", "-", names(expData))
    getNames <- function(x){ names(x) }
  } else {
    colnames(expData) <- gsub("\\.", "-", colnames(expData))
    getNames <- function(x) { colnames(x) }
  }
  rownames(phenodat) <- gsub("\\.", "-", rownames(phenodat))
  # ExpSamples <- bcIDR(getNames(expData))
  # duplicates <- as.logical(anyDuplicated(ExpSamples))
  # if (duplicates) {
  #   if (is(expData, "list")){
  #     expData <- expData[!duplicated(ExpSamples)]  
  #   } else {
  #     expData <- expData[!duplicated(ExpSamples),]
  #   }
  # }
  commonNames <- intersect(bcIDR(getNames(expData)), bcIDR(rownames(phenodat)))
  namesRight <- getNames(expData)[match(commonNames, bcIDR(getNames(expData)))]
  righttab <- bcRight(namesRight)
  clindup <- matrix(NA, nrow = length(commonNames))
  ## Rownames here use samples
  rownames(clindup) <- namesRight
  ## Match on patient identifiers
  clindup <- phenodat[match(commonNames, rownames(phenodat)),]
  clindup <- cbind(clindup, righttab[, -length(righttab)])
  ## Move patient IDs to a column (to allow any duplicates)
  clindup <- data.frame(patientids = rownames(clindup), clindup,
                         row.names = NULL, stringsAsFactors = FALSE)
  return(clindup)
}