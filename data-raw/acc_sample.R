library(RTCGAToolbox)

args <- formals(getFirehoseData)
setTRUE <- vapply(args,
    function(x) {
        is.logical(x) && !isTRUE(x)
    }, logical(1L)
)
args[setTRUE] <- TRUE
args$getUUIDs <- FALSE
args <- args[-length(args)]
args$dataset <- "ACC"

accmini <- do.call(getFirehoseData, args)

accmini@RNASeq2GeneNorm[[1]]@DataMatrix <- head(acc@RNASeq2GeneNorm[[1]]@DataMatrix)
accmini@miRNASeqGene <- head(acc@miRNASeqGene)
accmini@CNASNP <- head(acc@CNASNP)
accmini@CNVSNP <- head(acc@CNVSNP)
accmini@Methylation[[1]]@DataMatrix <- head(acc@Methylation[[1]]@DataMatrix)
accmini@RPPAArray[[1]]@DataMatrix <- head(acc@RPPAArray[[1]]@DataMatrix)
accmini@GISTIC@AllByGene <- head(acc@GISTIC@AllByGene)
accmini@GISTIC@ThresholdedByGene <- head(acc@GISTIC@ThresholdedByGene)
accmini@GISTIC@Peaks <- head(acc@GISTIC@Peaks)
accmini@Mutation <- head(acc@Mutation)

usethis::use_data(accmini, overwrite = TRUE)
