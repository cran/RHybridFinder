% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhc_check.R
\name{mhc_check}
\alias{mhc_check}
\title{mhc_check}
\usage{
mhc_check(netmhcpan_alleles)
}
\arguments{
\item{netmhcpan_alleles}{the netmhcpan alleles to be used for the netmhcpan call.}
}
\value{
\itemize{
          \item returns a custom error message if MHC/HLA allele(s) are not
          written correctly
          \item returns nothing if there are no issues. If HLA alleles are not
          written correctly}
}
\description{
this function only contains the alleles list, read by netMHCpan,
the list was retrieved by reading the file exported from netMHCpan, using the
following command line "netMHCpan -listMHC"
}
\details{
a custom error is printed in case the allele is not written correctly
}
\examples{
if (interactive()) {
 mhc_check("HLA-A02:01")
 mhc_check("HLA-A0201")
}
}
