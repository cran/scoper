# Project documentation and imports

#' The SCOPer package
#'
#' Provides a computational framework for unsupervised identification B cell
#' clones from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets.
#' This method is based on spectral clustering of the junction sequences of B cell
#' receptors (BCRs, Immunoglobulins) that share the same V gene, J gene and
#' junction length.
#'
#' @section Spectral Clustering for clOne Partitioning (SCOPer):
#'
#' \itemize{
#'   \item  \link{defineClonesScoper}:  Clustering sequences into clonal groups.
#'   \item  \link{analyzeClones}:       Summary statistics and visualization of the
#'                                      clonal clustering results.
#' }
#'
#' @name        scoper
#' @docType     package
#' @references
#' \enumerate{
#'   \item  Nouri N and Kleinstein SH (2018). A spectral clustering-based method for identifying clones
#'   from high-throughput B cell repertoire sequencing data. Bioinformatics, 34(13):i341-i349.
#' }
#'
#' @import      methods
#' @importFrom  ggplot2     ggplot aes_string
#'                          theme theme_bw element_text element_blank element_rect
#'                          ggtitle xlab ylab
#'                          scale_fill_manual
#'                          geom_polygon geom_histogram geom_hline geom_vline
#' @importFrom  dplyr       do n desc funs %>%
#'                          as_data_frame data_frame data_frame_
#'                          bind_cols bind_rows combine
#'                          filter filter_ select select_ arrange arrange_
#'                          group_by group_by_ ungroup
#'                          mutate mutate_ summarize summarize_
#'                          mutate_at summarize_at vars slice one_of
#' @importFrom  stringi     stri_split_fixed stri_length
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  alakazam    pairwiseDist checkColumns getDNAMatrix progressBar writeChangeoDb
#'                          groupGenes baseTheme
#' @importFrom  stats       density kmeans sd uniroot dnorm
#' @importFrom  iterators   icount
#' @importFrom  lazyeval    interp
NULL
