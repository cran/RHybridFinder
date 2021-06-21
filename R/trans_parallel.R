#' @title trans_parallel
#' @description non-exportable function, searches for peptide fragment combinations
#' existing in two parental proteins using parallel computing
#' @param nbCores number of cores to be used for parallel computing
#' @param fragments_not_cis_list the fragment combinations of non-cis peptides,
#' search_for_cis_peptides output
#' @param proteome_db proteome database
#' @return list of trans-spliced peptides
#' @details this function uses parallel computing in order to search for peptide
#' pair fragment existence in 2 proteins.(parallel version)
#' @noRd
#' @keywords internal
#' @seealso
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[foreach]{foreach}}
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach `%dopar%` foreach
#' @importFrom stats setNames


trans_parallel<- function(nbCores, fragments_not_cis_list, proteome_db){
  doParallel::registerDoParallel(nbCores)
  i<-NULL
  `%dopar%` <- foreach::`%dopar%`
  trans_spliced<- foreach::foreach(i=fragments_not_cis_list, .final=function(x) stats::setNames(x, names(fragments_not_cis_list))) %dopar% { find_trans_spliced_peptides_parallel(i, proteome_db)}
  doParallel::stopImplicitCluster()
  return (trans_spliced)
}
