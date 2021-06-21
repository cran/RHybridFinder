#' @title find_cis_spliced_peptides_parallel
#' @description this is a non-exportable function, the goal of this function is
#' to simply look for the fragments within proteins in the proteome database.
#' @param sequence_patt sequence patterns
#' @param proteome_db proteome database
#' @return list of the peptides that have a match for both fragments within the
#' same protein and that dont.
#' @details this is a non-exportable function, that takes in as input the sequence
#' patterns and the proteome database to search into and returns a list of all
#' all these patterns and their respective hits in the proteome, and blank if there
#' wasn't a hit.
#' @noRd
#' @keywords internal


find_cis_spliced_peptides_parallel<- function(sequence_patt, proteome_db){
  cis_spliced<- grep(sequence_patt, proteome_db, perl=TRUE)
  return(cis_spliced)
}
