#' @title find_trans_spliced_peptides
#' @description non-exportable functions that searches for fragments within 2
#' proteins.
#' @param x list of fragment combinations
#' @param proteome_db proteome database
#' @return a list of the trans-spliced peptides
#' @details this function searches the proteome for peptide fragments existing in
#' 2 proteins.
#' @noRd
#' @keywords internal

find_trans_spliced_peptides<- function(x, proteome_db){
  x<- setNames(strsplit(as.character(x), split="_"), names(x))
  s <- proteome_db
  trans_spliced_1<-mapply(function(y,z){
    if(length(grep(y, s, perl=TRUE))>0 & length(grep(z, s, perl=TRUE))>0){
      "trans"
    }else{
      "not trans"}
  }, y = sapply(x, "[", 1),z = sapply(x, "[", 2))
  trans_spliced_1<- setNames(trans_spliced_1, names(x))
  return(trans_spliced_1)
}
