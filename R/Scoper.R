# Project documentation and imports

#' The SCOPer package
#'
#' Provides a computational framework for B cell clones identification
#' from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets. 
#' Three models are included (identical, hierarchical, and spectral) 
#' which perform clustering among sequences of B cell receptors 
#' (BCRs, also referred to as Immunoglobulins, (Igs)) that 
#' share the same V gene, J gene and junction length.
#'
#' @section SCOPer:
#'
#' \itemize{
#'   \item  \link{defineClonesScoper}:  Clustering sequences into clonal groups.
#' }
#'
#' @name        scoper
#' @docType     package
#' @references
#' \enumerate{
#'   \item  Nouri N and Kleinstein SH (2018). A spectral clustering-based method for identifying clones
#'   from high-throughput B cell repertoire sequencing data. Bioinformatics, 34(13):i341-i349.
#'   \item Gupta NT, et al. (2017). Hierarchical clustering can identify B cell clones with high 
#'   confidence in Ig repertoire sequencing data. The Journal of Immunology, 198(6):2489-2499.
#' }
#'
#' @import      methods
#' @importFrom  ggplot2     ggplot aes_string 
#'                          theme theme_bw element_text element_blank element_rect
#'                          ggtitle xlab ylab coord_flip
#'                          scale_fill_manual scale_y_continuous geom_density
#'                          geom_polygon geom_histogram geom_hline geom_vline
#' @importFrom  dplyr       n %>% data_frame filter select arrange 
#'                          group_by ungroup group_indices
#'                          mutate summarize slice 
#' @importFrom  stringi     stri_split_fixed stri_length stri_count
#' @importFrom  stringr     str_count
#' @importFrom  seqinr      s2c
#' @importFrom  data.table  as.data.table .I
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  alakazam    pairwiseDist checkColumns getDNAMatrix padSeqEnds
#'                          progressBar groupGenes baseTheme translateDNA
#' @importFrom  shazam      consensusSequence
#' @importFrom  stats       density kmeans sd cor
#'                          as.dist hclust cutree
#' @importFrom  rlang       sym syms
#' @importFrom  Rcpp        evalCpp
#' @useDynLib   scoper, .registration=TRUE
NULL
