#' @title cis_parallel
#' @description non-exportable function, searches for peptide fragment combinations
#' within the same protein using parallel computing
#' @param nbCores number of cores to be used for parallel computing
#' @param isnot_Linear the fragment combinations of non-linear peptides,
#' search_for_linear_peptides output
#' @param proteome_db proteome database
#' @return list of cis-spliced peptides
#' @details this function uses parallel computing in order to search for cis-spliced
#' peptides
#' @keywords internal
#' @noRd
#' @seealso
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[foreach]{foreach}}
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach `%dopar%` foreach
#' @importFrom stats setNames

cis_parallel<- function(nbCores, isnot_Linear, proteome_db){
  # parallel computing only for finding cis-spliced peptides
  doParallel::registerDoParallel(nbCores)
  i<-NULL
  `%dopar%` <- foreach::`%dopar%`
  forward_cis_spliced<- foreach::foreach(i=isnot_Linear$splicePattern, .final=function(x) {stats::setNames(x, isnot_Linear$id )}) %dopar% { find_cis_spliced_peptides_parallel(i, proteome_db)}
  doParallel::stopImplicitCluster()
  return(forward_cis_spliced)
}
