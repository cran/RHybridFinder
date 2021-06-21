#' @title export_step2_results
#' @description this function allows to export the results generated from
#' step2_wo_netmhcpan.
#' @param step2_RHF_results_Exp1 the results generated from running
#' step2_wo_netMHCpan()
#' @param export_dir the export directory where you would like to have the .csv
#' file saved.
#' @return exports a folder containing 2 files
#' \enumerate{
#'           \item the peptide list to be entered in a netMHCpan-ready format,(.csv)
#'           \item the updated database search results which contain the categorizatiosn
#'            of the peptides found in common between the 2nd database search and
#'            the HybridFinder function (.csv file)}
#' @details Since netMHCpan is not compatible with Windows OS, the package offers
#' an alternative by outputting the input for netMHCpan and as well the database
#' results with their respective categorizations (cis, trans) established in step1.
#' @examples
#' \dontrun{
#'  export_step2_results(results_step2_Human_Liver_AUTD17, folder_Human_Liver_AUTD17)
#' }
#' @rdname export_step2_results
#' @export
#' @importFrom utils write.csv

export_step2_results <- function(step2_RHF_results_Exp1,export_dir){

  #retrieve dataframes
  peptide_list_nmhcp  <- step2_RHF_results_Exp1[[1]]
  rerun_db_updated <- step2_RHF_results_Exp1[[2]]

  #create naming convention
  obj_name<- deparse(substitute(step2_RHF_results_Exp1))
  unique_identifier<- grep(pattern="\\d+", Sys.time(), value=TRUE)
  unique_identifier<- gsub(pattern="\\:", "_", unique_identifier)

  #create directory
  unique_exportdir_name_Step2 <-  paste0(unique_identifier,"_", obj_name, "_step2_HF_results")
  output_dir<- file.path(export_dir,unique_exportdir_name_Step2 )
  dir.create(output_dir)

  # Write out final table containing hybrid peptides:
  utils::write.table(peptide_list_nmhcp, file.path(output_dir,'peptide_list_netMHCpan.pep'), row.names = FALSE,col.names = FALSE, quote = FALSE)
  utils::write.csv(rerun_db_updated,file = file.path(output_dir,'database_search_rerun.csv'),row.names = FALSE)

}
