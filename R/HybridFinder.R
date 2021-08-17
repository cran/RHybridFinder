#' @title HybridFinder
#' @description This function takes in three mandatory inputs: (1) all denovo
#' candidates (2) database search results and (3) the corresponding proteome
#' fasta file. The function's role is to extract high confidence de novo peptides
#' and to search for their existence in the proteome, whether the entire peptide
#' sequence or its pair fragments (in one or two proteins).
#' @param denovo_candidates dataframe containing all denovo candidate peptides
#' @param db_search dataframe containing the database search peptides
#' @param proteome_db path to the proteome FASTA file
#' @param customALCcutoff the default is calculated based on the median ALC of the
#' assigned spectrum groups (spectrum groups that match in the database search
#' results and in the denovo sequencing results) where also the peptide sequence
#' matches, Default: NULL
#' @param with_parallel for faster results, this function also utilizes
#' parallel computing (please read more on parallel computing in order
#' to be sure that your computer does support this), Default: TRUE
#' @param customCores custom amount of cores strictly higher than 5, Default: 6
#' @param export_files a boolean parameter for exporting the dataframes into
#' files in the next parameter for the output directory, Default: FALSE,
#' Default: FALSE
#' @param export_dir the output directory for the results files if
#' export_files=TRUE, Default: NULL, Default: NULL
#' @return The output is a list of 3 dataframes containing:
#' \enumerate{
#'            \item the HybridFinder output (dataframe) - the spectra that made
#'            it to the end with their respective columns (ALC, m/z, RT, Fraction,
#'            Scan) and a categorization column which denotes their potential splice
#'            type (-cis, -trans) or whether they are linear (the entire sequence
#'            was matched in proteins in the proteome database). Potential cis- &
#'            trans-spliced peptide are peptides whose fragments were matched with
#'            fragments within one protein, or two proteins, respectively.
#'            \item character vector containing potentially hybrid peptides (cis-
#'            and trans-)
#'            \item list containing the reference proteome and the "fake" proteins
#'            added at the end with a patterned naming convention (sp|denovo_HF_fake_protein)
#'            made up of the concatenated potential hybrid peptides.}
#' @details This function is based on the published algorithm by Faridi et al.
#' (2018) for the identification and categorization of hybrid peptides. The function
#' described here adopts a slightly modified version of the algorithm for
#' computational efficiency. The function starts by extracting unassigned denovo
#' spectra where the Average Local Confidence (assigned by PEAKS software), is
#' equivalent to the ALC cutoff which is based on the median of the assigned spectra
#' (between denovo and database search). The sequences of all peptides are searched
#' against the reference proteome. If there is a hit then, then, the peptide sequence
#' within a spectrum group considered as being linear and each spectrum group is
#' is then filtered so as to keep the highest ALC-ranking spectra. Then, the rest
#' of the spectra (spectra that did not contain any sequence that had an entire
#' match in the proteome database) then undergo a "cutting" procedure where each
#' sequence yields n-2 sequences (with n being the length of the peptide. That is
#' if the peptide contains 9 amino acids i.e NTYASPRFK, then the sequence is cut
#' into a combination of 7 sequences of 2 fragment pairs each i.e fragment 1: NTY
#' and fragment 2: ASPRFK, etc).These are then searched in the proteome for hits
#' of both peptide fragments within a same protein, spectra in which sequences
#' have fragment pairs that match within a same protein, these are considerent
#' to be potentially cis-spliced. Potentially cis-spliced spectrum groups are then
#' filtered based on the highest ranking ALC. Spectrum groups not considered to be
#' potentially cis-spliced are further checked for potential trans-splicing. The
#' peptide sequences are cut again in the same fashion, however, this time peptide
#' fragment pairs are searched for matches in two proteins. Peptide sequences whose
#' fragment pairs match in 2 proteins are considerend to be potentially trans-spliced.
#' The same filtering for the highest ranking ALC within each peptide spectrum group.
#' The remaining spectra that were neither assigned as linear nor potentially spliced
#' (neither cis- nor trans-) are then discarded. The result is a list of spectra
#' along with their categorizations (Linear, potentially cis- and potentially trans-)
#' Potentially cis- and trans-spliced peptides are then concatenated and then broken into
#' several "fake" proteins and added to the bottom of the reference proteome. The
#' point of this last step is to create a merged proteome (consisting of the reference
#' proteome and the hybrid proteome) which  would be used for a second database
#' search. After the second database search the checknetmhcpan function or the
#' step2_wo_netMHCpan function can be used in order to obtain the final list
#' of potentially spliced peptides.
#' Article: Faridi P, Li C, Ramarathinam SH, Vivian JP, Illing PT, Mifsud NA,
#' Ayala R, Song J, Gearing LJ, Hertzog PJ, Ternette N, Rossjohn J, Croft NP,
#' Purcell AW. A subset of HLA-I peptides are not genomically templated: Evidence
#' for cis- and trans-spliced peptide ligands. Sci Immunol. 2018 Oct 12;3(28):eaar3947.
#' <doi: 10.1126/sciimmunol.aar3947>. PMID: 30315122.
#' @examples
#' \dontrun{
#'  hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
#'  proteome, export = TRUE, output_dir)
#'  hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
#'  proteome)
#'  hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
#'  proteome, export = FALSE)
#' }
#' @seealso
#'  \code{\link[seqinr]{read.fasta}},\code{\link[seqinr]{s2c}}
#' @rdname HybridFinder
#' @export
#' @importFrom seqinr read.fasta s2c

HybridFinder<- function(denovo_candidates, db_search, proteome_db,
                        customALCcutoff=NULL, with_parallel= TRUE, customCores=6,
                        export_files = FALSE, export_dir = NULL){

  #use path for concatenating at the end
  proteome_path <- gsub(".*(/)","", proteome_db)


  #extract high quality and unique denovo peptides
  input_for_HF <- prepare_input_for_HF(denovo_candidates, db_search, customALCcutoff)

  #create an extra id column
  input_for_HF$extraid <- paste0(input_for_HF$denovo_id,"-+xx", input_for_HF$Peptide,
                                 "--xx", input_for_HF$ALC...., "_x",
                                 seq(1, nrow(input_for_HF)))

  #import proteome database
  proteome_db <- seqinr::read.fasta(proteome_db, as.string = TRUE, seqtype = "AA")

  #switch all "I" in the proteome to "L" since denovo does not distinguish between them
  proteome_db <- lapply(proteome_db, function(z) {gsub("I", "L", z)} )

  #search for linear peptides
  message('Step01: Search for linear peptides...')
  linear_peptides <- search_for_linear_peptides(input_for_HF, proteome_db,
                                                with_parallel, customCores)
  not_linear_peptides <- linear_peptides[[1]]
  final_df_linear_peptides <- linear_peptides[[2]]

  #search for cis-spliced peptides
  message('Step02: Search for cis-spliced peptides...')
  cis_spliced_peptides <- search_for_cis_spliced_peptides(not_linear_peptides ,
                                                          proteome_db,
                                                          with_parallel,
                                                          customCores)
  not_cis_peptides <- cis_spliced_peptides[[1]]
  final_df_cis_peptides <- cis_spliced_peptides[[2]]


  #search for trans-spliced peptides
  message('Step03: Search for trans-spliced peptides...')
  trans_spliced_peptides <- search_for_trans_spliced_peptides(not_cis_peptides,
                                                              proteome_db,
                                                              with_parallel,
                                                              customCores)
  final_df_trans_peptides <- trans_spliced_peptides

  #compile all peptides
  final_df_all_peptides <- rbind(final_df_linear_peptides, final_df_cis_peptides,
                                 final_df_trans_peptides )
  final_df_all_peptides <- final_df_all_peptides[
    which(!is.na(final_df_all_peptides$Peptide)),]


  #final_df_all_peptides<- final_df_all_peptides[,-2]
  metadata_denovo<- as.data.frame(matrix(unlist(strsplit(final_df_all_peptides$denovo_id, "-")),
                       nrow=nrow(final_df_all_peptides),
                       ncol=4, byrow=TRUE), stringsAsFactors = F)
  metadata_denovo[,3] <- as.numeric(metadata_denovo[,3])
  metadata_denovo[,4] <- as.numeric(metadata_denovo[,4])

  colnames(metadata_denovo) <- c("Fraction", "Scan", "m/z", "RT")

  #remove unnecessary columns
  cols_to_delete <- grep("Fragment|spliceType|denovo_id", colnames(final_df_all_peptides))
  final_df_all_peptides <- final_df_all_peptides[, -cols_to_delete]
  final_df_all_peptides<- cbind(metadata_denovo, final_df_all_peptides)

  #keep only 9-12 amino acids
  final_df_all_peptides <- final_df_all_peptides[nchar(final_df_all_peptides$Peptide) > 8,]
  final_df_all_peptides <- final_df_all_peptides[nchar(final_df_all_peptides$Peptide) < 13,]

  # only spliced and filter for 9-12 -mers
  only_spliced <- final_df_all_peptides[final_df_all_peptides$Type!= "Linear",
                                        c(grep("^Peptide$",
                                               colnames(final_df_all_peptides)),
                                          grep("^Length$",
                                               colnames(final_df_all_peptides)))]
  only_spliced <- only_spliced[!duplicated(only_spliced$Peptide),]

  #make fake proteins
  hybrid_proteome <- make_fake_HF_proteins(only_spliced)

  #concatenate proteome with hybrid proteins
  hybrid_concat <- list(proteome_db, hybrid_proteome)
  hybrid_concat <- do.call(c, hybrid_concat)

  #get only the peptide column
  only_spliced <- only_spliced$Peptide

  #remove the denovo_idpep column
  final_df_all_peptides$proteome_database_used <- proteome_path
  colnames(final_df_all_peptides)[grep("Type", colnames(final_df_all_peptides))] <- "Potential_spliceType"

  # return list of all peptides with categorization, the spliced peptides, the
  # cis-spliced peptides with details and the hybrid proteome
  results_list<- list(final_df_all_peptides, only_spliced, hybrid_concat)

  #export the results list into three files
  if (export_files == TRUE && !is.null(export_dir)) {
    if (dir.exists(export_dir)){
      export_HybridFinder_results(results_list, export_dir)
    }
  }
  #create the final output
  return (results_list)

}
