#' @title netmhcpan_call
#' @description non exportable function. netmhcpan_call formats the alleles and
#' runs netMHCpan.
#' @param hybrid_f output of peptide_rerun_cleanup, which is a  dataframe containing
#' the #' peptide sequences without modifications and of length 9-12 amino acids.
#' @param netmhcpan_alleles the alleles to be tested against with netMHCpan
#' @param path_for_tmp_files path to the temporary folder where lies the copied
#' netMHCpan file and where the output is saved.
#' @return does not return. Runs netMHCpan and outputs a results file in the
#' temporary folder.
#' @details this function simply formats the alleles and runs netMHCpan.
#' @noRd
#' @keywords internal
#' @importFrom utils write.table

netmhcpan_call<- function(hybrid_f, netmhcpan_alleles, path_for_tmp_files){

  utils::write.table(hybrid_f, 'peptides.pep', row.names = FALSE,
              col.names = FALSE, quote = FALSE)

  cat('running netMHCpan')

  unique_filename<- paste0(Sys.Date(),"_",deparse(substitute(hybrid_f)),
                           "_cnmhp_netmhcpan_output")

  #check if alleles are correctly written
  netmhcpan_alleles<- gsub("\\*", "", netmhcpan_alleles)
  mhc_check(netmhcpan_alleles)
  netmhcpan_alleles <- ifelse(length(netmhcpan_alleles) >1,noquote(paste(
    netmhcpan_alleles,collapse=",")), noquote(netmhcpan_alleles))
  #call netMHCpan
  system(paste0(file.path(path_for_tmp_files, "netMHCpan")," -a ",
                netmhcpan_alleles, " -p peptides.pep -inptype 1 -BA > ",
                unique_filename, ".pep"))

  nemhcpan_name<-file.path(paste0(unique_filename,".pep"))

  return(nemhcpan_name)
  }
