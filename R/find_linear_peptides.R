#' @title find_linear_peptides
#' @description this is a non-exportable function, the goal of this function is
#' to check for linear peptides among the high-ALC unassigned  spectra from the
#' denovo results.
#' @param input_for_HF the prepare_input_for_HF output
#' @param proteome_db proteome database
#' @return list of linear high-ALC denovo peptides and those that will be explored
#' further
#' @details this is a non-exportable function, that takes in as input the peptide
#' sequences and the proteome database to search into and returns a list of all
#' all these peptides and their respective hits in the proteome, and blank if there
#' wasn't a hit.
#' @noRd
#' @keywords internal
#' @importFrom seqinr getAnnot


find_linear_peptides<- function(input_for_HF, proteome_db){
  isLinear_peptides<- tapply(input_for_HF$Peptide,input_for_HF$extraid,function(x){
  seqinr::getAnnot(proteome_db[grep(x, proteome_db, fixed=TRUE)])})
return(isLinear_peptides)
}
