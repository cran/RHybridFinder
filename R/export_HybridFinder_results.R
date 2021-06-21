#' @title export_HybridFinder_results
#' @description this function allows to export the results list obtained in the
#' HybridFinder() function.
#' @param results_list the results list obtained with the HybridFinder() function.
#' @param export_dir the export directory
#' @return exports a folder containing three files
#' \enumerate{
#'            \item the HybridFinder output - the spectra that made it to the end
#'            with their respective columns (ALC, m/z, RT, Fraction, Scan) and a
#'            categorization column which denotes their potential splice type
#'            (-cis, -trans) or whether they are linear (the entire sequence was
#'            matched in proteins in the proteome database). Potential cis- &
#'            trans-spliced peptide are peptides whose fragments were matched with
#'            fragments within one protein, or two proteins, respectively.
#'            \item list of potential hybrid peptides (excluding the linear peptides)
#'            (.csv file)
#'            \item the merged proteome consisting of the reference proteome along
#'            with the hybrid proteome added at the end of the file with the sequence
#'            names following the pattern "sp|denovo_HF_fake_protein" along with
#'            a digit at the end (1,2,3,4,4,etc.) (.fasta file)}
#' @details In order to be able to have the HybridFinder() results list exported,
#' this function will come in handy. Please note that this function is also part
#' of the HybridFinder() function, therefore if you set export_files=TRUE and you
#' indicate the export directory in export_dir in the HybridFinder() function,
#' you would have the exact same outcome.
#' @examples
#' \dontrun{
#'  export_results(results_HybridFinder_Human_Liver_AUTD17,folder_Human_Liver_AUTD17)
#' }
#' @seealso
#'  \code{\link[seqinr]{write.fasta}}
#' @rdname export_results
#' @export
#' @importFrom seqinr write.fasta

export_HybridFinder_results <- function (results_list, export_dir) {

  #create a unique identifier
  unique_identifier<- grep(pattern="\\d+", Sys.time(), value=TRUE)
  unique_identifier<- gsub(pattern="\\:", "_", unique_identifier)
  obj_name<- deparse(substitute(results_list))

  #create output directory inside export_dir
  unique_exportdir_name_Step1<-  paste0(unique_identifier, "_", obj_name,"_", "HybrindFinder_Results")
  output_dir<- file.path(export_dir,unique_exportdir_name_Step1 )
  dir.create(output_dir)

  #write spliced peptides file
  write.table(as.data.frame(results_list[[2]], stringsAsFactors=F), file = file.path(output_dir,paste(obj_name,"spliced_peptides.csv", sep="_")), row.names = FALSE, quote = FALSE, col.names= FALSE)

  #write step05-final decision file
  hybrid_ori<- data.frame(lapply(results_list[[1]], as.character), stringsAsFactors=FALSE)
  write.csv(hybrid_ori, file = file.path(output_dir,"HF_output_spliced.csv"), row.names = FALSE, quote = FALSE, na="N/A")

  #create output filename for fasta
  proteome_filename <- ifelse(!is.null(obj_name),paste0("step1_HF_hybrid_proteome",deparse(substitute(results_list)), ".fasta"), "step1_HF_hybrid_proteome.fasta")

  #write concatenated fasta
  seqinr::write.fasta(sequences = results_list[[3]], names = names(results_list[[3]]), file.out =file.path(output_dir, proteome_filename))

}
