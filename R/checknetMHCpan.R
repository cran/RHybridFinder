#' @title checknetMHCpan
#' @description checknetMHCpan, utilizes the file from the second (PEAKS) run
#' and analyzes the data with netMHCpan in order to provide the peptide binding
#' affinity to different HLA/MHC alleles.
#' @param netmhcpan_directory the directory in which the netMHCpan file is.
#' @param netmhcpan_alleles vector of comma-separated alleles for which these peptides
#' should be analyzed (i.e HLA_alleles_Exp1<- c("HLA-A*02:01", "HLA-A*03:01",
#' "HLA-A24:02"))
#' @param peptide_rerun dataframe containing the results of the second run
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
#'    \item netMHCpan results pertaining to the binding affinity of all peptides
#' in the database search results (in long- and wide- format, with data tidying
#' in the wide format in order to compute the amount of HLA molecules to which
#' a peptide is strong/weak/non-binder binder)(dataframe)
#'     \item netMHCpan results pertaining to the binding affinity of the hybrid
#' peptides to the MHC molecules (in long- and wide- format, with data tidying
#' in the wide format in order to compute the amount of HLA molecules to which
#' a peptide is strong/weak/non-binder binder) (dataframe)
#'     \item the database search rerun with the categorizations already determined
#'     in step1 (HybridFinder Function)} (datafrane)
#' @details The ability to check the peptide binding affinity to the different
#' MHC/HLA molecules is essential for assessing the antigenicity of all peptides.
#' This function thus uses netMHCpan (Reynisson et al., 2020) for the generation
#' of binding affinty results.
#' @examples
#' \dontrun{
#'   results_checknetmhcpan_Exp1<- checknetMHCpan('/usr/local/bin', alleles,
#'   peptide_rerun, Exp1_HF_results[[1]])
#'   results_checknetmhcpan_Exp1 <- checknetMHCpan('/usr/local/bin', alleles,
#'   peptide_rerun, Exp1_HF_results_denovo_w_spliced)
#' }
#' @rdname checknetMHCpan
#' @export


checknetMHCpan <- function(netmhcpan_directory, netmhcpan_alleles, peptide_rerun,
                           HF_step1_output, export_files=FALSE,
                           export_dir=NULL){
  if (length(grep("ALC", colnames(peptide_rerun)))==0){

  }else{
    stop("Please provide the proper input")
  }

  user_prev_dir<- getwd()
  on.exit(setwd(user_prev_dir),add=TRUE)
  tryCatch(expr={
    if (export_files == FALSE && is.null(export_dir)) {
      export <- FALSE
      output_dir <- NULL
    }else if (export_files == TRUE && dir.exists(export_dir)==FALSE){
    stop("The directory you have provided is not a valid. Exiting.")
  }


  #create a temporary folder and temporary folder directory object
  path_for_tmp_files <- create_tmp_directory()
  on.exit(unlink(path_for_tmp_files), add=TRUE)
  setwd(path_for_tmp_files)
  Sys.sleep(1.5)
  cat("Temporary directory created:", path_for_tmp_files,"\n")
  Sys.sleep(1.5)


  # get default value in case peptide_col is not prrovided
  peptide_col<- grep("^Peptide$", colnames(peptide_rerun))

  #retrieve categorizations and prepare input for netMHCpan
  peptide_rerun_cleanup_list <- peptide_rerun_cleanup(peptide_rerun, HF_step1_output,
                                                      peptide_col)
  peptide_rerun <- peptide_rerun_cleanup_list[[1]]
  hybrid_f <- peptide_rerun_cleanup_list[[2]]

  #copy netMHCpan into tmp directory
  file.copy(from = grep( "netMHCpan$",list.files (path=netmhcpan_directory,
                                                  full.names=TRUE), value=TRUE),
            to = path_for_tmp_files)

  #calling netMHCpan
  netMHCpan_call <- netmhcpan_call(hybrid_f, netmhcpan_alleles, path_for_tmp_files)

  # getting the non-xls output of netMHCpan output and cleaning it
  netmhcpan_lists <- read_nmhcp_file(netMHCpan_call)
  df_netmhcpan_long<- netmhcpan_lists[[2]]
  df_netmhcpan_wide<- netmhcpan_lists[[1]]


  #get the potential splice type
  df_netmhcpan_long$Potential_spliceType<- peptide_rerun$Potential_spliceType[
    match(df_netmhcpan_long$Peptide, peptide_rerun$Peptide_no_mods)]

  df_netmhcpan_wide$Potential_spliceType<- peptide_rerun$Potential_spliceType[
    match(df_netmhcpan_wide$Peptide, peptide_rerun$Peptide_no_mods)]


  list_netMHCpan<- list(df_netmhcpan_long, df_netmhcpan_wide, peptide_rerun)
  #if the user would like to have the files exported
  if (export_files == TRUE && dir.exists(export_dir)){
    export_checknetmhcpan_results(list_netMHCpan,export_dir)
  }else{}

  #delete temporary directory
  unlink(path_for_tmp_files, recursive=TRUE)

  return(list_netMHCpan)},
  error=function(e){message(e)
    unlink(path_for_tmp_files, recursive = TRUE)
    },
  warning=function(w){message(w)
    unlink(path_for_tmp_files, recursive = TRUE)}
  )
}
