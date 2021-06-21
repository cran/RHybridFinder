#' @title find_trans_spliced_peptides_parallel
#' @description non-exportable functions that searches for fragments within 2
#' proteins.
#' @param x list of fragment combinations
#' @param proteome_db proteome database
#' @return a list of the trans-spliced peptides
#' @details this function searches the proteome for peptide fragments existing in
#' 2 proteins.
#' @noRd
#' @keywords internal

find_trans_spliced_peptides_parallel<- function(x, proteome_db){
  s <- proteome_db
  x<- setNames(strsplit(as.character(x), split="_"), names(x))

  if(length(grep(sapply(x, "[", 1), s, perl=TRUE))>0 & length(grep(
    sapply(x, "[", 2), s, perl=TRUE))>0){
      trans_spliced_1<-"trans"
    }else{
      trans_spliced_1<-"not trans"
    }
  return(trans_spliced_1)
}
