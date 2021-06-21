#' @title find_cis_start_positions
#' @description non-exportable function, goal is to find the positions of the
#' fragments within their parental protein
#' @param cis_spliced_df dataframe containing the cis spliced peptides along
#' with the patterns and indices
#' @param proteome_db proteome database
#' @return list containing all the fragment positions of the potentially
#' cis-spliced peptides.
#' @details this function finds the positions of the peptide fragments
#'  of potentially cis-spliced peptides within their parental protein in the
#'  proteome.
#' @noRd
#' @keywords internal


# function for finding fragment start positions
find_cis_start_positions<- function(cis_spliced_df, proteome_db){
  pattern <- cis_spliced_df$splicePattern
  index <- cis_spliced_df$indices
  s<- proteome_db
  positions<- mapply(function(x , y){
    attr(gregexpr(x, s[as.numeric(y)], perl=TRUE)[[1]],'capture.start')},x=pattern, y=index)

  return(positions)
  }
