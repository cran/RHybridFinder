#' @title mhc_check
#' @description this function only contains the alleles list, read by netMHCpan,
#' the list was retrieved by reading the file exported from netMHCpan, using the
#' following command line "netMHCpan -listMHC"
#' @param netmhcpan_alleles the netmhcpan alleles to be used for the netmhcpan call.
#' @details a custom error is printed in case the allele is not written correctly
#' @return
#' \itemize{
#'           \item returns a custom error message if MHC/HLA allele(s) are not
#'           written correctly
#'           \item returns nothing if there are no issues. If HLA alleles are not
#'           written correctly}
#' @examples
#' if (interactive()) {
#'  mhc_check("HLA-A02:01")
#'  mhc_check("HLA-A0201")
#' }
#' @rdname mhc_check
#' @export



mhc_check<- function(netmhcpan_alleles){

  #retrieve netMHCpan list of alleles
  netmhcpan_listofalleles<- netmhcpan_list_alleles$V1

  #re-format the input alleles
  if (length(grep("\\*", netmhcpan_alleles))==0){
    }else{
      netmhcpan_alleles<- gsub("\\*", "", netmhcpan_alleles)
    }

  #check if alleles are in netmhcpan list of alleles
  x <- sapply(strsplit(netmhcpan_alleles,split=","), function(y) {
    ifelse(y %in% netmhcpan_listofalleles, TRUE, FALSE)})
  if (all(x)) {
    #nothing returned if all good.
  }else{
    stop(" Please check the input alleles:", netmhcpan_alleles[x==FALSE])
  }
}
