#' Get resource links from inputs
#'
#' This function provides a reference to the resources downloaded from the
#' GDAC Firehose pipeline. Based on the input, the function returns a
#' URL location to the resource if there exists one.
#'
#' @inheritParams getFirehoseData
#'
#' @param data_date Either a runDate or analysisDate typically entered in
#' `getFirehoseData`
#'
#' @return A character URL to a dataset location
#' @examples
#'
#' getLinks("BRCA", CNASeq = TRUE)
#'
#' @export
getLinks <- function(dataset, data_date="20160128",
    RNASeqGene=FALSE, RNASeq2Gene=FALSE, clinical=FALSE, miRNASeqGene=FALSE,
    RNASeq2GeneNorm=FALSE,
    RNAseq2Norm = c(
      "normalized_counts", "RSEM_normalized_log2", "raw_counts", "scaled_estimate"
    ),
    CNASNP=FALSE, CNVSNP=FALSE, CNASeq=FALSE,
    CNACGH=FALSE, Methylation=FALSE, Mutation=FALSE, mRNAArray=FALSE,
    miRNAArray=FALSE, RPPAArray=FALSE, GISTIC=FALSE)
{
    fh_url <- "https://gdac.broadinstitute.org/runs/stddata__"

    if (GISTIC)
        fh_url <- "https://gdac.broadinstitute.org/runs/analyses__"

    fh_url <- paste0(fh_url, substr(data_date, 1, 4), "_",
        substr(data_date, 5, 6), "_", substr(data_date, 7, 8), "/data/")
    fh_url <- paste0(fh_url, dataset, "/", data_date, "/")

    doc <- .files_from_html_table(fh_url)

    paste0(
        fh_url,
    if (RNASeqGene)
        .getLinks("Level_3__gene_expression__data.Level_3", "*.Merge_rnaseq__.*._rnaseq__.*.tar[.]gz$", NULL, doc)
    else if (RNASeq2Gene)
        .getLinks("Level_3__RSEM_genes__data.Level_3", "*.Merge_rnaseqv2__.*._rnaseqv2__.*.tar[.]gz$", NULL, doc)
    else if (clinical)
        .getLinks(".Clinical_Pick_Tier1.Level_4", "*.tar[.]gz$", NULL, doc)
    else if (miRNASeqGene)
        .getLinks("Level_3__miR_gene_expression__data.Level_3", "[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$", dataset, doc)
    else if (RNASeq2GeneNorm) {
        RNAseq2Norm <- match.arg(RNAseq2Norm)
        if (!identical(RNAseq2Norm, "normalized_counts"))
            .getLinks("mRNAseq_Preprocess",".*\\.Level_3\\..*\\.tar\\.gz$",NULL,doc)
        else
            .getLinks("Level_3__RSEM_genes_normalized__data.Level_3", "*.Merge_rnaseqv2__.*._rnaseqv2__.*.tar[.]gz$", NULL, doc)
    } else if (CNASNP)
        .getLinks("Level_3__segmented_scna_hg19__seg.Level_3", "[.]Merge_snp__.*.__Level_3__segmented_scna_hg19__seg.Level_3.*.tar[.]gz$", dataset, doc)
    else if (CNVSNP)
        .getLinks("Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3", "[.]Merge_snp__.*.__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$", dataset, doc)
    else if (CNASeq)
        .getLinks("__Level_3__segmentation__seg.Level_3", "[.]Merge_cna__.*.dnaseq.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$", dataset, doc)
    else if (CNACGH)
        .getLinks("__Level_3__segmentation__seg.Level_3", "[.]Merge_cna__.*.cgh.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$", dataset, doc)
    else if (Methylation)
        .getLinks("__Level_3__within_bioassay_data_set_function__data.Level_3", "[.]Merge_methylation__.*.methylation.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$", dataset, doc)
    else if (Mutation)
        .getLinks("Mutation_Packager_Calls", "[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$", dataset, doc)
    else if (mRNAArray)
        c(
            .getLinks("Merge_transcriptome__agilentg4502a_07", "[.]Merge_transcriptome__agilentg4502a_.*.__Level_3__unc_lowess_normalization_gene_level__data.Level_3.*.tar[.]gz$", dataset, doc),
            .getLinks("Merge_transcriptome__ht_hg_u133a", "[.]Merge_transcriptome__ht_hg_u133a__.*.__Level_3__gene_rma__data.Level_3.*.tar[.]gz$", dataset, doc),
            .getLinks("Merge_exon__huex_1_0_st_v2", "[.]Merge_exon__huex_1_0_st_v2__.*.__Level_3__quantile_normalization_gene__data.Level_3.*.tar[.]gz$", dataset, doc)
        )
    else if (miRNAArray)
        .getLinks("h_mirna_8x15k", "[.]Merge_mirna__h_mirna_8x15k.*.data.Level_3.*.tar[.]gz$", dataset, doc)
    else if (RPPAArray)
        .getLinks("rppa_core", "[.]Merge_protein_exp.*.protein_normalization__data.Level_3.*.tar[.]gz$", dataset, doc)
    else if (GISTIC) {
        tag <- switch(dataset, SKCM = "TM", LAML = "TB", "TP")
        dset <- paste0(dataset, "-", tag)
        .getLinks("CopyNumber_Gistic2.Level_4", "[.]CopyNumber_Gistic2[.]Level_4.*.tar[.]gz$", dset, doc)
    }
    )
}

