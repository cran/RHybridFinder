#' @title search_for_linear_peptides
#' @description this is a non-exportable function. The goal of this function is
#' to find the denovo linear peptides.
#' @param input_for_HF the prepare_input_for_HF output
#' @param proteome_db proteome database
#' @param with_parallel whether parallel computing should be used, Default: TRUE
#' @return list of non-linear peptides and a table containing the linear peptides
#' formatted for the final rbind.
#' @details this function searches for the peptides in the proteome and returns
#' a table containing the linear peptides and a list containing the non-linear
#' peptides.
#' @noRd
#' @keywords internal
#' @examples
#' if (interactive()) {
#'  input_for_HF<- load_data_for_checking_examples ("test_input_for_HF")
#'  proteome_db <- load_data_for_checking_examples ("proteome_filtered")
#'  search_for_linear_peptides <- function(input_for_HF, proteome_db)
#' }
#' @importFrom stats setNames


#find linear peptides
search_for_linear_peptides <- function(input_for_HF, proteome_db, with_parallel, customCores){

  nbCores<- parallel::detectCores()

  if (with_parallel == TRUE & nbCores > 5){
    message('with parallel computing (cores)')
    nbCores<- customCores
    message(customCores)
    isLinear_peptides <- linear_parallel(nbCores, input_for_HF, proteome_db)
  } else if (with_parallel == TRUE & nbCores < 5){
    message('without parallel computing\n')
    isLinear_peptides<- find_linear_peptides(input_for_HF, proteome_db)
  } else if (with_parallel == FALSE){
    message('without parallel computing\n')
    isLinear_peptides<- find_linear_peptides(input_for_HF, proteome_db)
  }

  #retrieve the linear peptides
  linear_fd<- data.frame(Fragment=rep(names(isLinear_peptides),
                                      sapply(isLinear_peptides, length)),
                         fullannot=unlist(isLinear_peptides), stringsAsFactors = FALSE)

  #format the gene name
  linear_fd$GeneName <- gsub(".*GN=", "", linear_fd$fullannot)
  linear_fd$GeneName <-  gsub(" PE.*", "", linear_fd$GeneName)

  #retrieve denovo id
  linear_fd$denovo_id<- gsub("(\\-\\+xx).*","", linear_fd$Fragment)

  #create id column
  linear_fd$id <- gsub("(\\_x\\d+)$", "", linear_fd$Fragment)

  final_linear<- linear_fd

  #create ALC and peptide columns
  final_linear$ALC <- gsub(".*(\\-\\-xx)","", final_linear$id)
  final_linear$Peptide<- gsub(".*(\\-\\+xx)","", final_linear$id)
  final_linear$Peptide <- gsub("(\\-\\-xx).*","", final_linear$Peptide)
  final_linear$ALC<- as.numeric(final_linear$ALC)

  #get the linear peptides by keeping the maximum from each denovo id (spectrum)
  uniqulinearids<- unique(final_linear$denovo_id)
  final_linear$ALC <- as.numeric(final_linear$ALC)
  linear_spectrum<- by(final_linear, final_linear["denovo_id"], function(z) z[which(z$ALC[z$denovo_id %in% uniqulinearids] == max(z$ALC)),])
  linear_spectrum<- do.call(rbind, linear_spectrum)
  linear_spectrum <- linear_spectrum[!duplicated(linear_spectrum$id),]


  #find non linear for the next step
  not_linear <- data.frame(Names=names(isLinear_peptides[unlist(lapply(isLinear_peptides, function(x) length(x)<1))]),stringsAsFactors=FALSE)
  `%notin%` <- Negate(`%in%`)

  #create id, ALC, peptide, denovo id columns
  not_linear$id<- gsub("\\_x\\d+$","", not_linear$Names)
  not_linear$ALC <- gsub(".*(\\-\\-xx)","", not_linear$id)
  not_linear$Peptide<- gsub(".*(\\-\\+xx)","", not_linear$id)
  not_linear$Peptide <- gsub("(\\-\\-xx).*","", not_linear$Peptide)
  not_linear$denovo_id<- gsub("(\\-\\+xx).*","", not_linear$id)

  #keep the denovo ids that are not associated with linear peptides
  not_linear_ids<- not_linear$denovo_id %notin% linear_fd$denovo_id
  not_linear <- not_linear[not_linear_ids,]
  not_linear_list<- stats::setNames(as.list(not_linear$Peptide), not_linear$id)


  #start compiling the peptides - Linear peptides
  final_df_linear_peptides<- data.frame(Peptide=linear_spectrum$Peptide, Fragment=linear_spectrum$Peptide, denovo_id= linear_spectrum$denovo_id , Length=nchar(linear_spectrum$Peptide), Type="Linear", spliceType=NA, ALC = linear_spectrum$ALC, stringsAsFactors = FALSE)

  linear_peptides_output <- list(not_linear_list, final_df_linear_peptides)

  return(linear_peptides_output)
}
