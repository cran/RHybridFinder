#' @title get_detailed_cis_positions
#' @description this is a non-exportable function. The goal of this function is
#' to get a detailed result of the positions, gap, etc. of the cis-spliced peptides.
#' @param cis_spliced_fd table with the cis-spliced peptides
#' @param proteome_db proteome database
#' @return dataframe containing the cis-spliced peptides and their positions and gap
#' within the parental protein
#' @details this functions finds the start, end positions of the fragments of
#' cis-spliced peptides within their parental proteins
#' @noRd
#' @keywords internal


get_detailed_cis_positions<- function(cis_spliced_fd, proteome_db){

  #get the positions
  start_positions<- find_cis_start_positions(cis_spliced_df=cis_spliced_fd,
                                             proteome_db = proteome_db)

  if( is.list (start_positions)){

    new_start_positions<- convert_cis_positions_list_into_df(start_positions)
    new_start_positions$spliceType <- cis_spliced_fd$ForwardOrReverse

    #get the potential spliced candidate within a given protein
    amountof_matches_per_peptide_within_protein<- (ncol(new_start_positions)-2)/2

    #calculate the gap and keep peptides that have the lowest gap within the same protein
    for (i in seq(1,amountof_matches_per_peptide_within_protein)){
      p<- i*2
      k<- p-1
      new_start_positions[,paste0("difference",i)]<-
        #if gregexpr finds 0 for 2nd frag for forward, means they're not spliced
        ifelse(new_start_positions[,p]==0 & new_start_positions$spliceType=="forward",
               #if gregexpr finds 0 for 2nd frag for reverse, means 0 gap
               Inf, ifelse(new_start_positions[,p]==0 & new_start_positions$spliceType=="reverse",
                           #gregexpr provides the index of the previous letter for the 2nd frag
                           0, abs(new_start_positions[,p]- (new_start_positions[,k]+1))))
    }

    difference_cols<- grep("difference", colnames(new_start_positions))
    difference_df<- new_start_positions[,difference_cols]
    new_start_positions$min_difference<- apply(difference_df, 1, function(s)
    {
      colnames(difference_df)[which.min(s)]
    })
    new_start_positions$gap<-NULL

    for (i in seq(1,nrow(new_start_positions))){
      min_difference<- new_start_positions$min_difference
      new_start_positions$gap[i]<- new_start_positions[i,min_difference[i]]}
    new_start_positions$gap<- ifelse(new_start_positions$gap == Inf, NA, new_start_positions$gap)

    #create id before filtering
    new_start_positions$spliced_Peptide <- gsub("(\\?\\!.*\\?)", "_",new_start_positions$peptidePattern)
    new_start_positions$spliced_Peptide<- gsub("\\(|\\)", "",new_start_positions$spliced_Peptide)
    new_start_positions$spliced_fragDist_id<- paste(new_start_positions$spliced_Peptide, cis_spliced_fd$indices)

    new_start_positions<- new_start_positions[!is.na(new_start_positions$gap),]

    new_start_positions$min_difference_frag_column<- (as.numeric(gsub("difference", "",
                                                                      new_start_positions$min_difference))*2)-1
    new_start_positions$min_difference_frag_next<- new_start_positions$min_difference_frag_column+1
    for (i in seq(1,nrow(new_start_positions))){
      beststart<- new_start_positions$min_difference_frag_column
      new_start_positions$bestStart[i]<- new_start_positions[i,beststart[i]]}
    for (i in seq(1,nrow(new_start_positions))){
      bestend<- new_start_positions$min_difference_frag_next
      new_start_positions$bestEnd[i]<- new_start_positions[i,bestend[i]]}

    #create an id to pair back the dataframes
    cis_spliced_fd$spliced_fragDist_id<- paste(cis_spliced_fd$spliced_peptide, cis_spliced_fd$indices)

    #retrieve the gap
    cis_spliced_fd$bestStart<- new_start_positions$bestStart[match(cis_spliced_fd$spliced_fragDist_id, new_start_positions$spliced_fragDist_id)]
    cis_spliced_fd$bestEnd<- new_start_positions$bestEnd[match(cis_spliced_fd$spliced_fragDist_id, new_start_positions$spliced_fragDist_id)]
    cis_spliced_fd$gap<- new_start_positions$gap[match(cis_spliced_fd$spliced_fragDist_id, new_start_positions$spliced_fragDist_id)]

    #remove the id column
    cis_spliced_fd<- cis_spliced_fd[, -grep("spliced_fragDist_id", colnames(cis_spliced_fd))]

    new_start_positions<- cis_spliced_fd

  }

  else if (is.matrix(start_positions) & nrow(start_positions)==2 & ncol(start_positions) == nrow(cis_spliced_fd)){
    start_positions<- t(start_positions)
    new_start_positions<- cis_spliced_fd
    new_start_positions$bestStart <- start_positions[,1]
    new_start_positions$bestEnd <- start_positions[,2]
    new_start_positions$gap<- abs(new_start_positions$bestStart -
                                    new_start_positions$bestEnd)
    new_start_positions$gap<- ifelse(new_start_positions$gap == Inf, NA, new_start_positions$gap)
    new_start_positions<- new_start_positions[!is.na(new_start_positions$gap),]

  }

  #get the fragments' length
  new_start_positions$frag1Length <- sapply(new_start_positions$spliced_peptide, function(s) gregexpr("_", s, fixed=TRUE)[[1]]-1)
  new_start_positions$frag2Length  <- nchar(new_start_positions$spliced_peptide)-new_start_positions$frag1Length - 1

  #format the results
  almost_final_cis<- new_start_positions
  almost_final_cis$GeneName<- gsub(".*GN=", "", almost_final_cis$full_annot)
  almost_final_cis$GeneName<- gsub(" PE.*", "", almost_final_cis$GeneName)
  almost_final_cis$fullSeq<- gsub("\\_", "",almost_final_cis$spliced_peptide)
  frags<- matrix(unlist(strsplit(almost_final_cis$spliced_peptide, split="_")), ncol=2, byrow=TRUE)
  frags_df<- data.frame(fragment1=frags[,1], fragment2=frags[,2], stringsAsFactors = F)
  almost_final_cis<- cbind(almost_final_cis, frags_df)
  almost_final_cis$ReversetoForwardFullseq <- ifelse(almost_final_cis$ForwardOrReverse=="reverse", paste0(almost_final_cis$fragment2, almost_final_cis$fragment1), paste0(almost_final_cis$fragment1, almost_final_cis$fragment2))

  return (almost_final_cis)
}
