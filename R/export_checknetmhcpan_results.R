#' @title export_checknetmhcpan_results
#' @description this function allows to export the results generated from
#' checknetMHCpan()
#' @param list_checknetMHCpan_results the results generated from running
#' checknetMHCpan()
#' @param export_dir the export directory where the results .csv files should
#' be exported.
#' @return exports a folder containing three files
#' \enumerate{
#'            \item netMHCpan results in long format (the original output)(.csv file)
#'            \item netMHCpan results tidied (in wide format) so as to summarize
#'            the information per peptide (.tsv tab-separated file)
#'            \item the updated database search results which contain the categorizatiosn
#'            of the peptides found in common between the 2nd database search and
#'            the HybridFinder function (.csv file)}
#' @details In order to be able to have the checknetMHCpan() function results
#' exported, this function will come in handy. Please note that this function is
#' also part of the checknetMHCpan() function (if export_files is set to TRUE
#' and a valid export directory is indicated)
#' @examples
#' \dontrun{
#'  export_checknetmhcpan_results(results_checknetMHCpan_Human_Liver_AUTD17, folder_Human_Liver_AUTD17)
#' }
#' @rdname export_checknetmhcpan_results
#' @export
#' @importFrom utils write.csv

export_checknetmhcpan_results <- function(list_checknetMHCpan_results,export_dir){

  df_NetMHCpan_long  <- list_checknetMHCpan_results[[1]]
  df_NetMHCpan_wide <- list_checknetMHCpan_results[[2]]
  rerun_db_updated<- list_checknetMHCpan_results[[3]]
  obj_name<- deparse(substitute(list_checknetMHCpan_results))

  unique_identifier<- grep(pattern="\\d+", Sys.time(), value=TRUE)

  unique_identifier<- gsub(pattern="\\:", "_", unique_identifier)

  unique_exportdir_name_Step2 <-  paste0(unique_identifier,"_", obj_name, "_checknetMHCpan_results")

  output_dir<- file.path(export_dir,unique_exportdir_name_Step2 )

  dir.create(output_dir)

  # Write out final table containing hybrid peptides:
  utils::write.csv(df_NetMHCpan_long,file = file.path(output_dir,'checknetMHCpan_long_output.csv'),row.names = FALSE)
  utils::write.table(df_NetMHCpan_wide,file = file.path(output_dir,'checknetMHCpan_wide_output.tsv'),row.names = FALSE, quote=FALSE, sep="\t")
  utils::write.csv(rerun_db_updated,file = file.path(output_dir,'database_search_rerun.csv'),row.names = FALSE)

}
