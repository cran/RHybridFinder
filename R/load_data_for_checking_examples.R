#' @title load_data_for_checking_examples
#' @description function useful for loading data for the examples fields of many
#' funcitons.
#' @param object_needed name of the object (list/dataframe) needed
#' @return returns the object needed from all datasets for examples and for the
#' vignettes
#' @details function useful for loading data for the examples fields of many funcitons.
#' @keywords internal
#' @noRd

load_data_for_checking_examples<- function(object_needed){
  load('./R/sysdata.rda')
  object_needed <- mget(object_needed, ifnotfound = "not found")
  return(object_needed)
}
