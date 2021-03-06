% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_step2_results.R
\name{export_step2_results}
\alias{export_step2_results}
\title{export_step2_results}
\usage{
export_step2_results(step2_RHF_results_Exp1, export_dir)
}
\arguments{
\item{step2_RHF_results_Exp1}{the results generated from running
step2_wo_netMHCpan()}

\item{export_dir}{the export directory where you would like to have the .csv
file saved.}
}
\value{
exports a folder containing 2 files
\enumerate{
          \item the peptide list to be entered in a netMHCpan-ready format,(.csv)
          \item the updated database search results which contain the categorizatiosn
           of the peptides found in common between the 2nd database search and
           the HybridFinder function (.csv file)}
}
\description{
this function allows to export the results generated from
step2_wo_netmhcpan.
}
\details{
Since netMHCpan is not compatible with Windows OS, the package offers
an alternative by outputting the input for netMHCpan and as well the database
results with their respective categorizations (cis, trans) established in step1.
}
\examples{
\dontrun{
 export_step2_results(results_step2_Human_Liver_AUTD17, folder_Human_Liver_AUTD17)
}
}
