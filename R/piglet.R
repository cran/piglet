# piglet package documentation and import directives

#' The Program for Ig clusters (PIgLET) package
#'
#' PIgLET is a suite of computational tools that improves genotype inference and downstream AIRR-seq data analysis.
#' The package as two main tools. The first is Allele Clusters, this tool is designed to reduce the ambiguity within the IGHV alleles. The ambiguity
#' is caused by duplicated or similar alleles which are shared among different genes.
#' The second tool is an allele based genotype, that determined the presence of an allele based on
#' a threshold derived from a naive population.
#'
#' @section  Allele Similarity Cluster:
#' This section provides the functions that support the main tool of creating the allele similarity cluster form
#' an IGHV germline set.
#'
#' \itemize{
#'   \item  \link{inferAlleleClusters}:      The main function of the section to create the allele clusters based on a germline set.
#'   \item  \link{ighvDistance}:             Calculate the distance between IGHV aligned germline sequences.
#'   \item  \link{ighvClust}:                Hierarchical clustering of the distance matrix from `ighvDistance`.
#'   \item  \link{generateReferenceSet}:     Generate the allele clusters reference set.
#'   \item  \link{plotAlleleCluster}:        Plots the Hierarchical clustering.
#'   \item  \link{artificialFRW1Germline}:   Artificially create an IGHV reference set with framework1 (FWR1) primers.
#' }
#'
#' @section  Allele based genotype:
#' This section provides the functions to infer the IGHV genotype using
#' the allele based method and the allele clusters thresholds
#'
#' \itemize{
#'   \item  \link{inferGenotypeAllele}:      Infer the IGHV genotype using the allele based method.
#'   \item  \link{assignAlleleClusters}:     Renames the v allele calls based on the new allele clusters.
#'   \item  \link{germlineASC}:              Converts IGHV germline set to ASC germline set.
#'   \item  \link{recentAlleleClusters}:     Download the most recent version of the allele clusters table archive from zenodo.
#'   \item  \link{extractASCTable}:          Extracts the allele cluster table from the zenodo archive file.
#'   \item  \link{zenodoArchive}:            An R6 object to query the zenodo api.
#' }
#'
#' @name     piglet
#' @docType  package
#' @references
#' \enumerate{
#'   \item  ##
#'}
#'
#' @import   methods
#' @import   graphics
#' @import   grDevices
#' @import   utils
#' @import   dendextend
#' @import   ggplot2
#' @import   circlize
#' @import   jsonlite
#' @import   ComplexHeatmap
#' @importFrom  ggmsa            ggmsa facet_msa
#' @importFrom  R6               R6Class
#' @importFrom  dplyr            do n desc funs %>% distinct
#'                               as_data_frame data_frame
#'                               bind_cols bind_rows combine rowwise slice
#'                               filter select arrange
#'                               group_by ungroup
#'                               mutate mutate_ summarize summarize_
#'                               mutate_at summarize_at count_ count
#'                               rename rename_ transmute transmute_ pull ungroup row_number
#' @importFrom  data.table       := rbindlist data.table .N setDT CJ setorderv setkey .SD %chin%
#' @importFrom  stats            hclust as.dendrogram as.dist binom.test p.adjust setNames weighted.mean median
#' @importFrom  alakazam         getGene getAllele getFamily
#' @importFrom  rlang            .data
#' @importFrom  tigger           readIgFasta findUnmutatedCalls
#' @importFrom  Biostrings       DNAStringSet
#' @importFrom  DECIPHER         DistanceMatrix
#' @importFrom  RColorBrewer     brewer.pal.info brewer.pal
#' @importFrom  splitstackshape  cSplit
#' @importFrom  zen4R            download_zenodo
#' @importFrom  methods          setOldClass
#' @importFrom  splitstackshape  cSplit
NULL
