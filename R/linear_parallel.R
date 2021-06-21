#' @title linear_parallel
#' @description non-exportable function, searches for denovo linear peptides
#' @param input_for_HF the output of the "prepare_input_for_HF" function.
#' @param proteome_db proteome database
#' @return list of peptides that have matches in the proteome database and those
#' that dont.
#' @details this function searches for complete match of a peptide sequence in the
#' proteome database of the high confidence denovo peptides. (parallel version)
#' @noRd
#' @keywords internal
#' @seealso
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[foreach]{foreach}}
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach `%dopar%` foreach
#' @importFrom stats setNames

linear_parallel<- function(nbCores, input_for_HF, proteome_db){
  # parallel computing only for finding cis-spliced peptides
  doParallel::registerDoParallel(nbCores)
  i<-NULL
  `%dopar%` <- foreach::`%dopar%`
  linear_peptides<- foreach::foreach(i=input_for_HF$Peptide, .final=function(x) {stats::setNames(x, input_for_HF$extraid )}) %dopar% { find_linear_peptides_parallel(i, proteome_db)}
  doParallel::stopImplicitCluster()
  return(linear_peptides)
}
