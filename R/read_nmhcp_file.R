#' @title read_nmhcp_file
#' @description non exportable function. The goal of this function is to read the
#' generated netMHCpan output format it and analyse it
#' @param path_to_nmhcp_output the path to the netmhcpan file obtained with the
#' netmhcpan_call function
#' @return a list containing long- format dataframe of the netMHCpan
#' results (the regular netMHCpan output) and a wide- format where the data is
#' tidied up.
#' @details this function reads the netmhcpan file, formats it, creates 2 versions
#' a wide and a long format version. Also the function summarizes the results
#' in the strong/weak/none - binder and strong/weak/none binder_count columns. That
#' is, these columns represent the alleles to which each peptide would strongly/weakly
#' bind to as well as their respective count.
#' @noRd
#' @keywords internal
#' @importFrom utils read.table
#' @importFrom stats reshape


read_nmhcp_file<- function(path_to_nmhcp_output){
  df_NetMHCpan<- utils::read.table(text=(gsub("<= ", " ", readLines(path_to_nmhcp_output))), header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
  df_NetMHCpan<- df_NetMHCpan[-grep("HLA|H-2|H2", df_NetMHCpan$V1),]
  df_NetMHCpan<- df_NetMHCpan[-grep("[-]+", df_NetMHCpan$V1),]
  if (any(df_NetMHCpan$V1=="peptides")){df_NetMHCpan<- df_NetMHCpan[!df_NetMHCpan$V1=="peptides",]}else{}
  columnnames<- df_NetMHCpan[df_NetMHCpan$V3=="Peptide",][1,]
  df_NetMHCpan<- df_NetMHCpan[!df_NetMHCpan$V1=="Pos",]
  df_NetMHCpan<- df_NetMHCpan[!df_NetMHCpan$V1=="Protein",]
  if (any(df_NetMHCpan$V3=="")){df_NetMHCpan<- df_NetMHCpan[!df_NetMHCpan$V3=="",]}else{}
  if(any(grepl("of",df_NetMHCpan$V1))){df_NetMHCpan<- df_NetMHCpan[-grep("of",df_NetMHCpan$V1 ),]}else{}
  columnnames<- as.vector(unlist(columnnames[1,]))
  df_NetMHCpan[,ncol(df_NetMHCpan)]<-ifelse(df_NetMHCpan[,ncol(df_NetMHCpan)]=="SB", "Strong binder", ifelse(df_NetMHCpan[,ncol(df_NetMHCpan)]=="WB", "Weak binder", "Non binder"))
  if(columnnames[length(columnnames)]=="") {
    df_NetMHCpan<- df_NetMHCpan[, -grep("BindLevel$", colnames(df_NetMHCpan))]
    columnnames[length(columnnames)]<- "BindLevel"
  }else{}
  colnames(df_NetMHCpan)<- columnnames

  netMHCpanVersion<- ifelse(any(grep("%Rank_EL|%Rank_BA", colnames(df_NetMHCpan))),4.1, 4.0 )
  print(paste("you are using netMHCpan version", netMHCpanVersion))
  if(netMHCpanVersion== 4.1){
    Rank_HLA_EL<- grep("%Rank_EL", colnames(df_NetMHCpan))
    Rank_HLA_BA<- grep("%Rank_BA.", colnames(df_NetMHCpan))
  } else if(netMHCpanVersion == 4.0) {
    Rank_HLA<- grep("%Rank_|%Rank.|^%Rank$", colnames(df_NetMHCpan))
  }else{
  }


  # get binding alleles and the level per peptide
  mhc_column_index<- grep("^MHC$|^HLA$", colnames(df_NetMHCpan))
  df_NetMHCpan_o <- df_NetMHCpan
  df_NetMHCpan_o$strongBinder<- ifelse(df_NetMHCpan_o$BindLevel=="Strong binder", df_NetMHCpan_o[,mhc_column_index], "")
  df_NetMHCpan_o$weakBinder<- ifelse(df_NetMHCpan_o$BindLevel=="Weak binder", df_NetMHCpan_o[,mhc_column_index], "")
  df_NetMHCpan_o$noneBinder<- ifelse(df_NetMHCpan_o$BindLevel=="Non binder", df_NetMHCpan_o[,mhc_column_index], "")

  #create wide format
  if (netMHCpanVersion == 4.1) {
    df_NetMHCpan_wide<- data.frame(Peptide=df_NetMHCpan$Peptide, HLA=df_NetMHCpan[,mhc_column_index], Rank_BA=df_NetMHCpan[,Rank_HLA_BA] ,Rank_EL=df_NetMHCpan[,Rank_HLA_EL], BindLevel=df_NetMHCpan$BindLevel)
  } else if (netMHCpanVersion == 4.0){
    df_NetMHCpan_wide<- data.frame(Peptide=df_NetMHCpan$Peptide, HLA=df_NetMHCpan[, mhc_column_index], Rank=df_NetMHCpan[, Rank_HLA], BindLevel=df_NetMHCpan$BindLevel)
  }
  #df_NetMHCpan_wide<- data.frame(Peptide=df_NetMHCpan$Peptide, HLA=mhc_column_index, Rank=df_NetMHCpan$`%Rank`, BindLevel=df_NetMHCpan$BindLevel)
  df_NetMHCpan_wide$strongbinderAllele<- ifelse(grepl("Strong binder", df_NetMHCpan_wide$BindLevel),as.character(df_NetMHCpan_wide$HLA), "Not strong binder")
  df_NetMHCpan_wide$weakbinderAllele<- ifelse(grepl("Weak binder", df_NetMHCpan_wide$BindLevel),as.character(df_NetMHCpan_wide$HLA), NA)
  df_NetMHCpan_wide$binderAllele<- ifelse(df_NetMHCpan_wide$BindLevel != "Non binder",as.character(df_NetMHCpan_wide$HLA), NA)

  #create wide version
  pre_df_NetMHCpan_wide_1 <- data.frame(HLA=df_NetMHCpan_wide[,mhc_column_index], Peptide=df_NetMHCpan_wide$Peptide, strongbinderAllele=df_NetMHCpan_wide$strongbinderAllele)
  df_NetMHCpan_wide_1 <- stats::reshape(pre_df_NetMHCpan_wide_1, idvar="Peptide", timevar= grep("^MHC$|^HLA$", colnames(pre_df_NetMHCpan_wide_1),value=TRUE), direction="wide")

  dfm<- data.frame(strongBinder=tapply(df_NetMHCpan_o$strongBinder,df_NetMHCpan_o$Peptide, paste, collapse=","),
                   weakBinder=tapply(df_NetMHCpan_o$weakBinder,df_NetMHCpan_o$Peptide, paste, collapse=","),
                   noneBinder=tapply(df_NetMHCpan_o$noneBinder,df_NetMHCpan_o$Peptide, paste, collapse=","))
  dfm$Peptide<- row.names(dfm)
  dfm<- data.frame(dfm,row.names = NULL)

  dfm[,1]<-  gsub("\\,,+", ",", dfm[,1])
  dfm[,1]<-  gsub("^\\,$|^\\,|\\,$| ", "", dfm[,1])
  dfm[,2]<-  gsub("\\,,+", ",", dfm[,2])
  dfm[,2]<-  gsub("^\\,$|^\\,|\\,$| ", "", dfm[,2])
  dfm[,3]<-  gsub("\\,,+", ",", dfm[,3])
  dfm[,3]<-  gsub("^\\,$|^\\,|\\,$| ", "", dfm[,3])


  df_NetMHCpan <-  merge(x=df_NetMHCpan, y=dfm, by.x="Peptide", by.y="Peptide")
  mhc_column_index<- grep("^MHC$|^HLA$", colnames(df_NetMHCpan))
  #reshape to get HLA ranks in columns
  if (netMHCpanVersion == 4.1){
    pre_df_NetMHCpan_wide<- df_NetMHCpan[,c(1, mhc_column_index, Rank_HLA_EL)]
    df<-stats::reshape(pre_df_NetMHCpan_wide, idvar="Peptide", timevar=grep("^MHC$|^HLA$", colnames(pre_df_NetMHCpan_wide), value=TRUE), direction="wide")
  }else if (netMHCpanVersion == 4.0){
    pre_df_NetMHCpan_wide<- df_NetMHCpan[,c(1, mhc_column_index, Rank_HLA:ncol(df_NetMHCpan))]
    df<-stats::reshape(pre_df_NetMHCpan_wide, idvar="Peptide", timevar=grep("^MHC$|^HLA$", colnames(pre_df_NetMHCpan_wide), value=TRUE), direction="wide")
  }
  get_col_index<- grep("\\%Rank.|\\%Rank_|Peptide", colnames(df))
  df<- df[,get_col_index]

  df_NetMHCpan<- df_NetMHCpan[,-grep("Pos|Core|Of|Gp|Gl|Ip|Il|Icore|Identity|%Rank|Score|Aff(nM)", colnames(df_NetMHCpan))]
  df_NetMHCpan<- merge(x=df_NetMHCpan, y=df, by.x="Peptide", by.y="Peptide")

  # clean up the dataframe to have wide format and remove unnecessary columns
  df_NetMHCpan<- df_NetMHCpan[!duplicated(df_NetMHCpan$Peptide),]
  cols_to_delete<- grep("\\%Rank$|^HLA$|^MHC$|BindLevel|Aff", colnames(df_NetMHCpan))
  df_NetMHCpan<- df_NetMHCpan[,-cols_to_delete]

  # create columns for counting the alleles for which a peptide is a strong/weak/none binder
  df_NetMHCpan$strongBinder_count <-  lengths(strsplit(df_NetMHCpan$strongBinder, split=","))
  df_NetMHCpan$weakBinder_count<- lengths(strsplit(df_NetMHCpan$weakBinder, split=","))
  df_NetMHCpan$noneBinder_count<- lengths(strsplit(df_NetMHCpan$noneBinder, split=","))

  df_NetMHCpan_list<- list(df_NetMHCpan, df_NetMHCpan_o)

  return(df_NetMHCpan_list)
}
