#' @title step2_wo_netMHCpan
#' @description This function helps retrieve the categorizations for the peptides
#' from step 1 and apply them to those that are matched in the second database
#' search.
#' @param peptide_rerun dataframe containing the results of the second database
#' PEAKS search.
#' @param HF_step1_output the HybridFinder output containing the potential
#' splicing categorizations obtained with the HybridFinder function (HybridFinder)
#' based on the matching of fragment pairs of peptides in 1 or 2 proteins. This
#' parameter can be provided either by loading the .csv exported file, or if the
#' results #' object still is in the global environment (i.e results_HF_Exp1),
#' then it can be accessed by simply writing "results_HF_Exp1[[1]]".
#' @param export_files a boolean parameter for exporting the dataframes into
#' files in the next parameter for the output directory, Default: FALSE
#' @param export_dir export_dir the output directory for the results files
#' if export_files=TRUE, Default: NULL
#' @return
#' \enumerate{
#'     \item the input file for the web version of netMHCpan (dataframe)
#'     \item the database search rerun with the categorizations already determined
#'     in the previous step. (character vector)}
#' @details In special cases where the PC runs on windows OS, since it would only
#' be possible to use the web version of netMHCpan, this function returns the
#' peptide input file for the webversion of netMHCpan. Also, this function outputs
#' the database search rerun results with their categorizations (into potentially
#' cis and potentially trans) obtained from the first step (HybridFinder).
#' @examples
#' if (interactive()) {
#' data(package="RHybridFinder", "db_rerun_Human_Liver_AUTD17")
#' results_checknetmhcpan_Human_Liver_AUTD17<- step2_wo_netMHCpan(db_rerun_Human_Liver_AUTD17,
#'     results_HybridFinder_Human_Liver_AUTD17[[1]])
#' }
#' @rdname step2_wo_netMHCpan
#' @export


step2_wo_netMHCpan <- function(peptide_rerun, HF_step1_output,
                               export_files=FALSE, export_dir=NULL){

  if (length(grep("ALC", colnames(peptide_rerun)))!=0){
    stop("Please provide the proper input")
  }else{
  }
  # get default value in case peptide_col is not prrovided
  peptide_col<- grep("^Peptide$", colnames(peptide_rerun))

  #loading the peptides after user validates list
  peptide_rerun_cleanup_list <- peptide_rerun_cleanup(peptide_rerun, HF_step1_output,
                                                      peptide_col)
  peptide_rerun <- peptide_rerun_cleanup_list[[1]]
  hybrid_f <- peptide_rerun_cleanup_list[[2]]

  list_step2_results<- list(hybrid_f, peptide_rerun)
  #if the user would like to have the files exported
  if (export_files == TRUE && dir.exists(export_dir)){
    export_step2_results(list_step2_results,export_dir)
  }else{}

  return(list_step2_results)
}
