#' @title convert_cis_positions_list_into_df
#' @description convert list containing the cis positions into a df, non-exportable
#' function. Special case where traditional list to df conversion leads to
#' errors in orders of the columns
#' @param cis_positions_list list containing fragment 1, 2 positions
#' @return dataframe containing the fragments 1,2 positions
#' @details this function convert list containing the cis positions into a df,
#' non-exportable function. Special case to solve issue where
#' traditional list to df conversion leads to errors in orders of the columns.
#' @noRd
#' @keywords internal


convert_cis_positions_list_into_df<- function(cis_positions_list){
  cis_positions_df<- as.data.frame(matrix(nrow=length(cis_positions_list), ncol=max(lengths(cis_positions_list))))
  cis_positions_df$peptidePattern <- names(cis_positions_list)

   for (i in seq(1,length(cis_positions_list))){
         pr<- length(cis_positions_list[[i]])
         if (pr ==2){
               for (k in seq(1,pr)){
                 cis_positions_df[i,k] <- cis_positions_list[[i]][k]
                 }
         }else{
           cis_positions_df[i,1]<- cis_positions_list[[i]][1]
                 for (k in seq(2,pr)){
                       if(k%%2==0){
                             if(k ==pr ){
                               cis_positions_df[i,k]=cis_positions_list[[i]][k]
                               }else{
                                 cis_positions_df[i,k]=cis_positions_list[[i]][k+1]
                             }}else{
                               cis_positions_df[i,k]=cis_positions_list[[i]][k-1]
                                  }

                     }
             }
     }
  return (cis_positions_df)
}
