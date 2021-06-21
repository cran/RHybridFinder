#' @title create_spliced_peptides
#' @description this is a non-exportable function. The goal of this function is
#' to simply create all fragment combinations possible from non-linear peptides.
#' @param i list of non-linear peptides to be 'cut' in order to create all possible
#' combinations.
#' @return list containing all combinations possible of non-linear peptides where,
#' with consideration for both splice types: forward and reverse.
#' @details this function creates a list of all combinations possible of fragments
#' from a given list of non-linear peptides. It creates both forward and reverse
#' peptide fragment combinations where each fragment is no less than 2 amino acids.
#' @noRd
#' @keywords internal


create_spliced_peptides<- function(i){
  spliced<- lapply(i, function(pep){
    lapply(seq(2, nchar(pep)-2),function(y){
      a <- substr(pep, start=0, stop=y)
      b <- substr(pep, start=y+1, stop=nchar(pep))
      splicing<- c(list(peptide= pep, Frag1=a, Frag2=b, spliceType="forward",
                        spliced_peptide=paste(a,b,sep="_")), list(peptide=pep,
                                                                  Frag1=b, Frag2 =a, spliceType="reverse", spliced_peptide=paste(b,a,sep="_")))
      return (splicing)
    })
  })
  x <- names(unlist(spliced, recursive = T))
  x <- x[-grep("\\.Frag1|\\.Frag2|\\.spliceType|\\.spliced_peptide", x)]
  counts<-seq(1, length(x))
  x <- gsub("\\.peptide",".x" , x)
  x<- paste0(x, counts)
  spliced_df<- data.frame(matrix(unlist(spliced), ncol=5, byrow=TRUE), stringsAsFactors = FALSE,row.names= x)
  colnames(spliced_df)<- c("Peptide", "Fragment1", "Fragment2", "spliceType", "spliced_peptide")
  return(spliced_df)
}
