#' @title prepare_input_for_HF
#' @description this is a non-exportable function. The goal of this function is
#' to extract unique denovo peptides (not matched in the database search), and that
#' are high confidence (they have a high ALC, determined based on the
#' the median ALC for linear assigned spectra where the peptide sequences in denovo
#' matches that in database search.
#' @param denovo_candidates the all denovo candidates .csv imported file results
#' of the denovo sequencing.
#' @param db_search the all database search .csv imported file results of the
#' database search results
#' @return a list containing
#' \enumerate{
#'         \item a dataframe containing all the high quality & unique (which have
#'         not been matched with database search peptide) denovo peptides along
#'         with their denovo id.
#'          \item the denovo dataframe in order to merge it at the end with the
#'          potential spliced results.}
#' @details this function extracts unique and high confidence denovo peptides from
#' the denovo sequencing results and database search results, and prepares a table
#' to be analyzed further on
#' @noRd
#' @importFrom stats sd median

prepare_input_for_HF <- function(denovo_candidates, db_search){

  file_mz_denovoCol<- grep("m[[:punct:]]z",colnames(denovo_candidates))
  file_mz_dbCol <- grep("m[[:punct:]]z",colnames(db_search))
  file_alc_denovoCol <- grep("ALC", colnames(denovo_candidates))
  file_alc_denovoName <- grep("ALC", colnames(denovo_candidates), value=TRUE)
  if (length(file_alc_denovoCol)<1){
    stop("Please make sure that you have the right input. N.B: The denovo results dataframe should be the first input")
  }else{
    file_alc_denovoCol <- grep("ALC", colnames(denovo_candidates))
    file_alc_denovoName <- grep("ALC", colnames(denovo_candidates), value=TRUE)
  }

  # creating specific ID (fraction-Scan-mz-RT)
  denovo_candidates$denovo_id <- paste(denovo_candidates$Fraction, denovo_candidates$Scan, denovo_candidates[,file_mz_denovoCol], denovo_candidates$RT, sep = "-")
  db_search$db_id <- paste(db_search$Fraction, db_search$Scan, db_search[,file_mz_dbCol], db_search$RT, sep = "-")

  # Replacing all the I to L in the DB_search peptide column and pasting results new column (Peptide_ItoL)
  db_search$Peptide_ItoL <- gsub("I", "L", db_search$Peptide)

  # copy PeptideItoL from db_search into denovo_candidates if ID is similar in both denovo_candidates and db_search
  denovo_candidates$PeptideItoL_match <- db_search$Peptide_ItoL [match(denovo_candidates$denovo_id, db_search$db_id)]

  # match denovo peptide with db_search PeptideItoL and selecting all positive matches for ALC calculation in the next step
  denovo_candidates$Peptide_match_test <- ifelse(denovo_candidates$Peptide == denovo_candidates$PeptideItoL_match, 1, 0)
  positive_ALC <- denovo_candidates$Peptide_match_test == "1"

  denovo_candidates[, file_alc_denovoCol] <- as.numeric(denovo_candidates[, file_alc_denovoCol])

  # Calculating new ALC [median]
  Positive_ALC_median <- stats::median(denovo_candidates[positive_ALC,file_alc_denovoCol], na.rm = TRUE)
  #Positive_ALC_sd <- stats::sd(denovo_candidates[positive_ALC,file_alc_denovoCol], na.rm = TRUE)
  #cutoff_ALC<- Positive_ALC_median - Positive_ALC_sd
  cutoff_ALC <- Positive_ALC_median
  new_ALC <- ceiling(cutoff_ALC)

  # sorting out denovo only peptides based on new calculated ALC.
  denovo_only_ALC <- denovo_candidates [((denovo_candidates[,file_alc_denovoCol] >= new_ALC) & is.na(denovo_candidates$Peptide_match_test)), ]

  # Select only required columns
  denovo_only_ALC_2 <- denovo_only_ALC[c("denovo_id", "Peptide", file_alc_denovoName)]

  # removing mods
  denovo_only_ALC_2$Peptide <- gsub(" *\\(.*?\\) *", "", denovo_only_ALC_2$Peptide)

  # Remove duplicate
  denovo_only_ALC_2_nodup <- denovo_only_ALC_2[!duplicated(denovo_only_ALC_2), ]

  input_for_HF<- denovo_only_ALC_2_nodup

  return(input_for_HF)

}
