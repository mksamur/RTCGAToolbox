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
  namesRight <- getNames(expData)[bcIDR(getNames(expData)) %in% commonNames]
  righttab <- bcRight(namesRight)
  clindup <- matrix(NA, nrow = length(commonNames))
  ## Rownames here use samples
  rownames(clindup) <- namesRight
  ## Match on patient identifiers
  clindup <- phenodat[rownames(phenodat) %in% commonNames,]
  ##
  righttab <- righttab[match(rownames(clindup), bcIDR(righttab$patientids)),]
  clindup <- cbind(clindup, righttab[, -length(righttab)])
  clindup <- data.frame(patientids = rownames(clindup), clindup,
                         row.names = NULL)
  return(clindup)
}