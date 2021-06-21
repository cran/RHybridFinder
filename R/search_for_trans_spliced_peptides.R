#' @title search_for_trans_spliced_peptides
#' @description non-exportable function. Finds the peptides whose fragment pairs
#' match in 2 proteins
#' @param not_cis_peptides list of peptides that are not cis, out from the
#' search_for_cis_spliced_peptides
#' @param proteome_db proteome database
#' @param with_parallel whether parallel computing should be used, Default: TRUE
#' @return table containing peptides whose fragment pairs match in two proteins
#' thereby classified as potentially trans-spliced
#' @details this function determines which peptides have fragment pairs that can
#' match in 2 proteins. Those that don't are then considered as having NBE(no
#' Biological Explanation) and are discarded.
#' @noRd
#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom stats setNames


search_for_trans_spliced_peptides<- function(not_cis_peptides, proteome_db,  with_parallel){

   fragments_not_cis<- create_spliced_non_cis_peptides(not_cis_peptides)

   fragments_not_cis_list<-  stats::setNames(as.list(fragments_not_cis$spliced_peptide_tbc),
                                      rownames(fragments_not_cis))

   #check for the amount of cores
   nbCores<- parallel::detectCores()

   #look for the peptides in the proteome with or without parallel computing
   if (with_parallel == TRUE & nbCores > 5){
      message('with parallel computing\n')
     trans_spliced <- trans_parallel(nbCores, fragments_not_cis_list, proteome_db)
     trans_spliced_fd<- data.frame(Fragment=rep(names(trans_spliced),
                                                sapply(trans_spliced, length)),
                                   transorNot=unlist(trans_spliced),
                                   stringsAsFactors = FALSE)
   }else if (with_parallel == TRUE & nbCores < 5){
      message("without parallel computing\n")
      trans_spliced <- data.frame(find_trans_spliced_peptides(fragments_not_cis_list,
                                                              proteome_db))
      colnames(trans_spliced)[1] <- c("transorNot")
      trans_spliced_fd <- data.frame(Fragment= rownames(trans_spliced),
                                     transorNot = trans_spliced[,1],
                                     stringsAsFactors = FALSE)
   }else if (with_parallel == FALSE){
      message("without parallel computing\n")
      trans_spliced <- data.frame(find_trans_spliced_peptides(fragments_not_cis_list,
                                                              proteome_db))
      colnames(trans_spliced)[1] <- c("transorNot")
      trans_spliced_fd <- data.frame(Fragment= rownames(trans_spliced),
                                     transorNot = trans_spliced[,1],
                                     stringsAsFactors = FALSE)
   }

   #retrieve the trans spliced and unpack the columns
   trans_spliced_fd$Fragment<- gsub("(\\.x\\d+)", "", trans_spliced_fd$Fragment)
   trans_spliced_fd$Peptide<- gsub(".*(\\-\\+xx)" ,"", trans_spliced_fd$Fragment, perl=TRUE)
   trans_spliced_fd$ALC<- gsub(".*(\\-\\-xx)" ,"", trans_spliced_fd$Peptide, perl=TRUE)
   trans_spliced_fd$Peptide<- gsub("(\\-\\-xx).*" ,"", trans_spliced_fd$Peptide, perl=TRUE)
   trans_spliced_fd$ALC<- gsub("(\\-\\+\\*\\&x).*" ,"", trans_spliced_fd$ALC, perl=TRUE)
   trans_spliced_fd$denovo_id<- gsub("(\\-\\+xx).*" ,"", trans_spliced_fd$Fragment, perl=TRUE)
   trans_spliced_fd_transonly<-  trans_spliced_fd[trans_spliced_fd$transorNot=="trans",]
   trans_spliced_fd_transonly <- trans_spliced_fd_transonly[!duplicated(trans_spliced_fd_transonly$Fragment),]

   #keep only the maximum ALC from each denovo id (spectrum)
   uniquetransids<- unique(trans_spliced_fd_transonly$denovo_id)
   trans_spliced_fd_transonly$ALC <- as.numeric(trans_spliced_fd_transonly$ALC)
   trans_spliced_fd_transonly_interim <- by(trans_spliced_fd_transonly,
                                            trans_spliced_fd_transonly["denovo_id"],
                                            function(z) z[which(z$ALC[z$denovo_id %in% uniquetransids] == max(z$ALC)),])
   trans_spliced_fd_transonly_interim <- do.call(rbind, trans_spliced_fd_transonly_interim)

   #start compiling peptides - trans-spliced peptides
   trans_spliced_fd_final<- trans_spliced_fd_transonly_interim
   final_df_trans_peptides<- data.frame(Peptide=trans_spliced_fd_final$Peptide,
                                        Fragment=trans_spliced_fd_final$Peptide,
                                        denovo_id= trans_spliced_fd_final$denovo_id ,
                                        Length=nchar(trans_spliced_fd_final$Peptide),
                                        Type="trans", spliceType=NA, ALC = trans_spliced_fd_final$ALC,
                                        stringsAsFactors = FALSE)


   return(final_df_trans_peptides)

   }
