## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load package, eval=FALSE, echo=TRUE--------------------------------------
#  library(RHybridFinder)

## ----get example data, eval=FALSE, echo=TRUE----------------------------------
#  #retrieve the denovo sequencing results for the example data
#  data(package="RHybridFinder", denovo_Human_Liver_AUTD17)
#  
#  #retrieve the database search results for the example data
#  data(package="RHybridFinder", db_Human_Liver_AUTD17)

## ----load inputs, eval=FALSE, echo=TRUE---------------------------------------
#  folder_Human_Liver_AUTD17 <- file.path("./data/Human_Liver_AUTD17")
#  denovo_Human_Liver_AUTD17 <- read.csv(file.path(folder_Human_Liver_AUTD17, "first_run","all de novo candidates.csv"), sep=",", head=TRUE,stringsAsFactors = FALSE)
#  db_Human_Liver_AUTD17 <- read.csv(file.path(folder_Human_Liver_AUTD17, "first_run","DB search psm.csv"), sep=",", head=TRUE,stringsAsFactors = FALSE)
#  proteome_Human_Liver_AUTD17<- file.path(folder_Human_Liver_AUTD17, "uniprot-proteome-human_UP000005640-reviewed_validated.fasta")

## ----loading inputs, eval=TRUE, echo=FALSE------------------------------------
load("../R/sysdata.rda")

## ----run HybridFinder, eval=FALSE, echo=TRUE----------------------------------
#  
#  results_HybridFinder_Human_Liver_AUTD17<- HybridFinder(denovo_candidates =  denovo_Human_Liver_AUTD17, db_search =  db_Human_Liver_AUTD17, proteome_db = proteome_Human_Liver_AUTD17, with_parallel = FALSE, export_files = TRUE, export_dir=folder_Human_Liver_AUTD17)

## ----show first output, eval=TRUE, echo=TRUE, paged.print=TRUE----------------
#display HybridFinder(HF) step1 output
print(head(results_HybridFinder_Human_Liver_AUTD17[[1]]))


## ----show second output, eval=TRUE, echo=TRUE, paged.print=TRUE---------------
#display list of candidate hybrid peptides
print(head(results_HybridFinder_Human_Liver_AUTD17[[2]]))

## ----show third output, eval=TRUE, echo=TRUE, paged.print=TRUE----------------
#display the merged proteome
print(tail(results_HybridFinder_Human_Liver_AUTD17[[3]]))

## ----loading inputs for checknetMHCpan, eval=FALSE, echo=TRUE-----------------
#  netmhcpan_dir<- '/usr/bin/'
#  
#  alleles_Human_liver_AUTD17<- c("HLA-A*03:01", "HLA-A*24:02", "HLA-B*35:03", "HLA-B*45:01", "HLA-C*04:01", "HLA-C*16:01")
#  
#  db_rerun_Human_liver_AUTD17 <- read.csv(file.path(folder_Human_Liver_AUTD17, "second_run","DB search psm.csv"), sep=",", head=TRUE,stringsAsFactors = FALSE)
#  
#  HF_output_Human_liver_AUTD17<- results_HybridFinder_Human_Liver_AUTD17[[1]]
#  

## ----run checknetMHCpan, eval=FALSE, echo=TRUE--------------------------------
#  
#  results_checknetMHCpan_Human_Liver_AUTD17<- checknetMHCpan(netmhcpan_directory = netmhcpan_dir, netmhcpan_alleles = alleles_Human_liver_AUTD17, peptide_rerun = db_rerun_Human_liver_AUTD17, HF_step1_output = HF_output_Human_liver_AUTD17, export_files = TRUE, export_dir=folder_Human_Liver_AUTD17)

## ----show netmhcpan_long_output, eval=TRUE, echo=TRUE, paged.print=TRUE-------
#display netmhcpan output(long version)
print(head(results_checknetMHCpan_Human_Liver_AUTD17[[1]]))


## ----show netmhcpan_wide_output, eval=TRUE, echo=TRUE, paged.print=TRUE-------
#display netmhcpan output tidied version (wide)
print(head(results_checknetMHCpan_Human_Liver_AUTD17[[2]]))


## ----show database_updated_output, eval=TRUE, echo=TRUE, paged.print=TRUE-----
#display the updated database search results with the categorizations from step1
print(head(results_checknetMHCpan_Human_Liver_AUTD17[[3]]))


## ----loading inputs for step2_wo_netMHCpan, eval=FALSE, echo=TRUE-------------
#  db_rerun_Human_liver_AUTD17 <- read.csv(file.path(folder_Human_Liver_AUTD17, "second_run","DB search psm.csv"), sep=",", head=TRUE,stringsAsFactors = FALSE)
#  
#  HF_output_Human_liver_AUTD17<- results_HybridFinder_Human_Liver_AUTD17[[1]]
#  

## ----running step2_wo_netMHCpan, eval=FALSE, echo=TRUE------------------------
#  
#  results_step2_Human_Liver_AUTD17<- step2_wo_netMHCpan(peptide_rerun = db_rerun_Human_liver_AUTD17, HF_step1_output = HF_output_Human_liver_AUTD17, export_files = TRUE, export_dir=folder_Human_Liver_AUTD17)

## ----show netmhcpan-ready input, eval=TRUE, echo=TRUE, paged.print=TRUE-------
#display the netmhcpan-ready input / list of all peptides 9-12 aa, without 
#modifications
print(head(results_step2_Human_Liver_AUTD17[[1]]))


## ----show database_updated_output(2), eval=TRUE, echo=TRUE, paged.print=TRUE----
#display the updated database search results table with the categorizations from 
#step1
print(head(results_step2_Human_Liver_AUTD17[[2]]))


