#' @title search_for_cis_spliced_peptides
#' @description non-exportable function. Searches for peptide fragments that co-exist
#' in the same protein. when there are is more than one match for the same peptide
#' fragment pair within the same protein, the one with the shortest intervening
#' distance is the one that is kept.
#' @param not_linear_list output from search_for_linear_peptides
#' @param proteome_db proteome database
#' @param with_parallel whether parallel computing should be used, Default: TRUE
#' @return a list containing the hybrid peptides whose pair fragments did not
#' match in one protein and a table containing hybrid peptides that have 1 or more
#' matches of the pair fragment within a given protein and therefore are categorized
#' as potentially cis-spliced.
#' @details this function Searches for peptide fragments that co-exist in the
#' same protein. when there are is than one match for the same peptide fragment
#' pair within the same protein, the one with the shortest intervening distance
#' is the one that is kept.
#' @noRd
#' @keywords internal
#' @importFrom seqinr getAnnot
#' @importFrom parallel detectCores
#' @importFrom stats setNames


search_for_cis_spliced_peptides <-  function(not_linear_list, proteome_db, with_parallel){

  #create all possible fragment combinations
  isnot_Linear <- create_spliced_peptides(not_linear_list)
  isnot_Linear$id<- gsub("\\.x\\d+", "", rownames(isnot_Linear))
  isnot_Linear$splicePattern <-  paste0("(",isnot_Linear$Fragment1,
                                        "((?!",isnot_Linear$Fragment1, ").)*?",
                                        isnot_Linear$Fragment2,")")
  isnot_Linear$id <- paste0(isnot_Linear$id, "-+*&x", isnot_Linear$spliceType,
                            "-,&x", isnot_Linear$splicePattern)

  #check for the amount of cores
  nbCores<- parallel::detectCores()

  #look for the peptides in the proteome with or without parallel computing
  if (with_parallel == TRUE & nbCores > 5){
    message('with parallel computing\n')
    cis_spliced_peptides <- cis_parallel(nbCores, isnot_Linear, proteome_db)
  } else if (with_parallel == TRUE & nbCores < 5){
    message('without parallel computing\n')
    isnot_Linear_list<- stats::setNames(as.list(isnot_Linear$splicePattern), isnot_Linear$id)
    cis_spliced_peptides <- find_cis_spliced_peptides (isnot_Linear_list, proteome_db)
  } else if (with_parallel == FALSE){
    message('without parallel computing\n')
    isnot_Linear_list<- stats::setNames(as.list(isnot_Linear$splicePattern), isnot_Linear$id)
    cis_spliced_peptides <- find_cis_spliced_peptides (isnot_Linear_list, proteome_db)
  }

  #retrieve the cis spliced and unpack the id column
  cis_spliced_fd<- data.frame(Fragment=rep(names(cis_spliced_peptides),
                                           sapply(cis_spliced_peptides, length)),
                              indices=unlist(cis_spliced_peptides),
                              stringsAsFactors = FALSE)
  cis_spliced_fd$ForwardOrReverse<- gsub(".*(\\-\\+\\*\\&x)","", cis_spliced_fd$Fragment)
  cis_spliced_fd$ForwardOrReverse<- gsub("(\\-\\,\\&x).*","", cis_spliced_fd$ForwardOrReverse)
  cis_spliced_fd$Peptide<- gsub(".*(\\-\\+xx)" ,"", cis_spliced_fd$Fragment, perl=TRUE)
  cis_spliced_fd$ALC<- gsub(".*(\\-\\-xx)" ,"", cis_spliced_fd$Peptide, perl=TRUE)
  cis_spliced_fd$Peptide<- gsub("(\\-\\-xx).*" ,"", cis_spliced_fd$Peptide, perl=TRUE)
  cis_spliced_fd$ALC<- gsub("(\\-\\+\\*\\&x).*" ,"", cis_spliced_fd$ALC, perl=TRUE)
  cis_spliced_fd$denovo_id<- gsub("(\\-\\+xx).*" ,"", cis_spliced_fd$Fragment, perl=TRUE)
  cis_spliced_fd$splicePattern<- gsub(".*(\\-\\,\\&x)", "", cis_spliced_fd$Fragment)
  cis_spliced_fd$spliced_peptide<- gsub("(\\(\\(\\?.*\\?)", "_", cis_spliced_fd$splicePattern)
  cis_spliced_fd$spliced_peptide<- gsub("\\(|\\)", "", cis_spliced_fd$spliced_peptide)


  #get the origin proteins and the accession
  cis_spliced_fd$full_annot <-  seqinr::getAnnot(proteome_db[as.numeric(cis_spliced_fd$indices)])
  cis_spliced_fd$accessionNumber<- names(proteome_db[as.numeric(cis_spliced_fd$indices)])


  #prepare data for finding the start position of each pattern within the protein
  cis_spliced_fd$frag1 <- gsub("_.*", "", cis_spliced_fd$spliced_peptide)
  cis_spliced_fd$frag2 <- gsub(".*_", "", cis_spliced_fd$spliced_peptide)
  cis_spliced_fd$id<- paste0(cis_spliced_fd$denovo_id,"++",
                             cis_spliced_fd$indices,",,,+",
                             cis_spliced_fd$ForwardOrReverse,"[,.$-+",
                             cis_spliced_fd$spliced_peptide,"-+]",
                             cis_spliced_fd$ALC,"{++",
                             cis_spliced_fd$full_annot,"++}",
                             "+++", cis_spliced_fd$accessionNumber )

  #get the positions, intervening distance(gao) of the cis-spliced peptides
  almost_final_cis<- get_detailed_cis_positions(cis_spliced_fd, proteome_db)

  cis_spectrum_interim <- almost_final_cis

  #retrieve the cis-spliced peptides by keeping the maximum ALCs from each denovo
  #id (spectrum)
  cis_spectrum_interim$duplid <- paste(cis_spectrum_interim$denovo_id,
                                       cis_spectrum_interim$ReversetoForwardFullseq,
                                       cis_spectrum_interim$GeneName)
  cis_spectrum_interim <- cis_spectrum_interim[!duplicated(cis_spectrum_interim$duplid),]
  cis_spectrum_interim$duplid <- paste(cis_spectrum_interim$denovo_id,
                                       cis_spectrum_interim$ReversetoForwardFullseq)
  cis_spectrum_interim <- cis_spectrum_interim[!duplicated(cis_spectrum_interim$duplid),]
  uniquecisids<- unique(cis_spectrum_interim$denovo_id)
  cis_spectrum_interim$ALC <- as.numeric(cis_spectrum_interim$ALC)
  cis_spectrum_interim2<- by(cis_spectrum_interim, cis_spectrum_interim["denovo_id"],
                             function(z) z[which(z$ALC[z$denovo_id %in% uniquecisids] == max(z$ALC)),])
  cis_spectrum_interim2<- do.call(rbind, cis_spectrum_interim2)

  #start compiling peptides - cis spliced peptides
  final_cis<- cis_spectrum_interim2
  final_df_cis_peptides<- data.frame(Peptide=final_cis$ReversetoForwardFullseq,
                                     Fragment=final_cis$Peptide,
                                     denovo_id= final_cis$denovo_id ,
                                     Length=nchar(final_cis$ReversetoForwardFullseq), Type="cis",
                                     spliceType=final_cis$ForwardOrReverse,
                                     ALC= final_cis$ALC,
                                     stringsAsFactors = FALSE)



  #find non cis-spliced peptides
  not_cis <- data.frame(Names=names(cis_spliced_peptides[unlist(lapply(
    cis_spliced_peptides, function(x) length(x)<1))]),stringsAsFactors=FALSE)
  not_cis$id<- gsub("\\_x\\d+$","", not_cis$Names)
  not_cis$ALC <- gsub(".*(\\-\\-xx)","", not_cis$id)
  not_cis$Peptide<- gsub(".*(\\-\\+xx)","", not_cis$id)
  not_cis$Peptide <- gsub("(\\-\\-xx).*","", not_cis$Peptide)
  not_cis$denovo_id<- gsub("(\\-\\+xx).*","", not_cis$id)
  not_cis$ALC <- gsub("(\\-\\+\\*\\&x).*","", not_cis$ALC)
  not_cis_ids<- not_cis$denovo_id %in% cis_spectrum_interim2$denovo_id
  not_cis <- not_cis[!not_cis_ids,]
  not_cis$id <- gsub("(\\-\\+\\*\\&x).*", "", not_cis$id)
  not_cis<- not_cis[!duplicated(not_cis$id),]
  not_cis_list<- stats::setNames(as.list(not_cis$Peptide), not_cis$id)
  not_cis_list <- not_cis_list[which(!is.na(not_cis_list))]

  cis_peptides_output<- list(not_cis_list, final_df_cis_peptides)

  return(cis_peptides_output)

}
