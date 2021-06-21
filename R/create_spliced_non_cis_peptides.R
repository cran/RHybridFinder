#' @title create_spliced_non_cis_peptides
#' @description non exportable function. this function creates fragment
#' combinations in order to search for trans-spliced peptides
#' @param i non cis- spliced peptides, output of search_for_cis_spliced_peptides
#' function
#' @return list of fragment combinations
#' @details this function creates fragment combinations from non-potentially cis
#' spliced peptides list in order to search for potential trans-spliced peptides.
#' @noRd
#' @keywords internal


create_spliced_non_cis_peptides<- function(i){
  spliced<- lapply(i, function(pep){
    lapply(seq(2, nchar(pep)-2),function(y){
      a <- substr(pep, start=0, stop=y)
      b <- substr(pep, start=y+1, stop=nchar(pep))
      splicing<- list(peptide= pep, Frag1=a, Frag2=b, spliceType="not cis", spliced_peptide_tbc =paste(a, b, sep="_") )
      return (splicing)
    })
  })
  x <- names(unlist(spliced, recursive = T))
  x <- x[-grep("\\.Frag1|\\.Frag2|\\.spliceType|\\.spliced_peptide_tbc", x)]
  counts<-seq(1, length(x))
  x <- gsub("\\.peptide",".x" , x)
  x<- paste0(x, counts)
  spliced_df<- data.frame(matrix(unlist(spliced), ncol=5, byrow=TRUE), stringsAsFactors = FALSE,row.names= x)
  colnames(spliced_df)<- c("Peptide", "Fragment1", "Fragment2", "spliceType", "spliced_peptide_tbc")
  return(spliced_df)
}
