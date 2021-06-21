#' @title make_fake_HF_proteins
#' @description This is a non-exportable function. The function is inpired by the
#' makefakeproteins from Faridi et al(2018)'s article for creating a fasta file
#' by concatenating the spliced peptides found.
#' @param only_spliced 8-12 mer spliced peptides
#' @return a list of concatenated peptides
#' @details this function concatenates the candidate hybrid peptides in a
#' write.fasta-ready format. It concatenates maximum 40 peptides per protein.The
#' objective of this function is to generate a hybrid proteome that would be
#' in the second step of HybridFinder which consists in running a second database
#' search using the merged proteome (reference proteome put in + the hybrid
#' proteome)
#' @noRd
#' @keywords internal
#' @importFrom stats setNames

make_fake_HF_proteins<-  function(only_spliced){

  #split the spliced peptides dataframe
  if (nrow(only_spliced)>40){
    create_proteins_40aa <- split(only_spliced[,1], cut(seq(nrow(only_spliced)), ceiling(nrow(only_spliced)/40)))
    #concatenate per split
    merge_peptides <- sapply(create_proteins_40aa, paste0, collapse="")
  }else{
    create_proteins_40aa <- only_spliced[,1]
    merge_peptides <- paste0(create_proteins_40aa, collapse="")
  }

  #create fasta names
  denovo_fake_protein_names <- paste0('sp|denovo_HF_fake_protein', seq(1, length(merge_peptides)))

  #attach the names to the list
  merge_peptides<- stats::setNames(merge_peptides, denovo_fake_protein_names)

  hybrid_proteome <-as.list(merge_peptides)
  return (hybrid_proteome)
}
