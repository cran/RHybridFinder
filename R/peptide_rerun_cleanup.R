#' @title peptide_rerun_cleanup
#' @description This function takes the dataframe containing the second database
#' search results.The function's goal is to retrieve categorizations from step1,
#' remove potential bias, prepare the input for netMHCpan. This is a non-exportable
#' function.
#' @param peptide_rerun the dataframe containing the second database search
#' results
#' @param HF_step1_output the HybridFinder output containing the potential
#' splicing categorizations obtained with the HybridFinder function (HybridFinder)
#' based on the matching of fragment pairs of peptides in 1 or 2 proteins.
#' @param peptide_col the number of the column containing the peptide sequences.
#' @param lower_limit_mer the lower limit in terms of peptide sequence length
#' or in other words, the minimum amount of amino acids per sequence, Default: 9
#' @param upper_limit_mer the upper limit in terms of peptide sequence length
#' or in other words, the maximum amount of amino acids per sequence, Default: 12
#' @return list containing:
#' \enumerate{
#'    \item the database search results updated with the appropriate categorizations
#'    for the different peptide sequences based on their identification in step1.
#'    Also, from the resulting dataframe are excluded all peptides that matched
#'    solely with the fake proteins but were not declared as cis/trans in step1,
#'    thereby removing any potential bias resulting from the addition of the hybrid
#'    proteome to the reference proteome.
#'    \item peptide sequences without modification and filtered for a length that
#'    can be used in netMHCpan (9-12 aa)}
#' @details This function updates the database search results with appropriate
#' categorizations from step 1 and removes the modifications from the peptide
#' sequences of the second database search results and filters peptides in order
#' to obtain peptides between 9 and 12 amino acids.
#' @noRd
#' @keywords internal


peptide_rerun_cleanup<- function(peptide_rerun, HF_step1_output , peptide_col, lower_limit_mer=9, upper_limit_mer=12) {

  # For universal accessibility, selecting the peptide column
  colnames(peptide_rerun)[peptide_col]<- "Peptide"

  #merge the denovo dataframe from HybridFinder (step1) with database rerun results
  HF_step1_output <- as.data.frame(HF_step1_output)
  HF_step1_output$Potential_spliceType <- as.factor(HF_step1_output$Potential_spliceType)
  HF_step1_output$Potential_spliceType<- factor(HF_step1_output$Potential_spliceType, levels=c("cis", "trans", "Linear"))
  HF_step1_output<- HF_step1_output[order(ordered(HF_step1_output$Peptide), HF_step1_output$Potential_spliceType),]

  peptide_rerun$Peptide_no_mods <- gsub(" *\\(.*?\\) *", "", peptide_rerun$Peptide)
  peptide_rerun$Potential_spliceType<- HF_step1_output$Potential_spliceType[match(peptide_rerun$Peptide_no_mods,
                                                                                  HF_step1_output$Peptide)]

  peptide_rerun$Potential_spliceType[
    is.na(peptide_rerun$Potential_spliceType)]<- "Linear"


  #remove peptides that are found only in fake proteins and were not in step1 output
  peptide_rerun$new_acc<- gsub("\\|denovo_HF_fake_protein\\d+||denovo_HF_fake_protein\\d+\\:|\\:|denovo_HF_fake_protein\\d+",
                               "", peptide_rerun$Accession)
  peptide_rerun$new_acc<- as.factor(peptide_rerun$new_acc)
  `%notin%` <- Negate(`%in%`)
  fake_biased_indx<- which(peptide_rerun$Potential_spliceType =="Linear" &
                             peptide_rerun$Peptide_no_mods %notin% HF_step1_output$Peptide &
                             peptide_rerun$new_acc=="")

  peptide_rerun<- peptide_rerun[-fake_biased_indx,-ncol(peptide_rerun)]

  #select useful info
  hybrid_f<- peptide_rerun$Peptide_no_mods
  hybrid_f<- hybrid_f[nchar(hybrid_f)>=lower_limit_mer]
  hybrid_f<- hybrid_f[nchar(hybrid_f)<=upper_limit_mer]
  hybrid_f<- unique(hybrid_f)

  peptide_rerun_cleanup_list <- list(peptide_rerun, hybrid_f)

  return (peptide_rerun_cleanup_list)
}
