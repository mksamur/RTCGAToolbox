#' RTCGAToolbox: A New Tool for Exporting TCGA Firehose Data
#'
#' Managing data from large-scale projects (such as The Cancer Genome Atlas (TCGA) for further analysis is an important and time consuming step for research projects. Several efforts, such as the Firehose project, make TCGA pre-processed data publicly available via web services and data portals, but this information must be managed, downloaded and prepared for subsequent steps. We have developed an open source and extensible R based data client for pre-processed data from the Firehose, and demonstrate its use with sample case studies. Results show that our RTCGAToolbox can facilitate data management for researchers interested in working with TCGA data. The RTCGAToolbox can also be integrated with other analysis pipelines for further data processing.
#'
#' The main function you're likely to need from `RTCGAToolbox` is
#' [getFirehoseData]. Otherwise refer to the vignettes to see
#' how to use the `RTCGAToolbox`
#'
#' @author Mehmet Kemal Samur
#' @name RTCGAToolbox
"_PACKAGE"

#' Gene coordinates for circle plot.
#'
#' A dataset containing the gene coordinates
#' The variables are as follows:
#'
#' * GeneSymbol: Gene symbols
#' * Chromosome: Chromosome name
#' * Strand: Gene strand on chromosome
#' * Start: Gene location on chromosome
#' * End: Gene location on chromosome
#'
#' @format A data frame with 28454 rows and 5 variables
#' @name hg19.ucsc.gene.locations
NULL

