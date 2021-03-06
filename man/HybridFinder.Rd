% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HybridFinder.R
\name{HybridFinder}
\alias{HybridFinder}
\title{HybridFinder}
\usage{
HybridFinder(
  denovo_candidates,
  db_search,
  proteome_db,
  customALCcutoff = NULL,
  with_parallel = TRUE,
  customCores = 6,
  export_files = FALSE,
  export_dir = NULL
)
}
\arguments{
\item{denovo_candidates}{dataframe containing all denovo candidate peptides}

\item{db_search}{dataframe containing the database search peptides}

\item{proteome_db}{path to the proteome FASTA file}

\item{customALCcutoff}{the default is calculated based on the median ALC of the
assigned spectrum groups (spectrum groups that match in the database search
results and in the denovo sequencing results) where also the peptide sequence
matches, Default: NULL}

\item{with_parallel}{for faster results, this function also utilizes
parallel computing (please read more on parallel computing in order
to be sure that your computer does support this), Default: TRUE}

\item{customCores}{custom amount of cores strictly higher than 5, Default: 6}

\item{export_files}{a boolean parameter for exporting the dataframes into
files in the next parameter for the output directory, Default: FALSE,
Default: FALSE}

\item{export_dir}{the output directory for the results files if
export_files=TRUE, Default: NULL, Default: NULL}
}
\value{
The output is a list of 3 dataframes containing:
\enumerate{
           \item the HybridFinder output (dataframe) - the spectra that made
           it to the end with their respective columns (ALC, m/z, RT, Fraction,
           Scan) and a categorization column which denotes their potential splice
           type (-cis, -trans) or whether they are linear (the entire sequence
           was matched in proteins in the proteome database). Potential cis- &
           trans-spliced peptide are peptides whose fragments were matched with
           fragments within one protein, or two proteins, respectively.
           \item character vector containing potentially hybrid peptides (cis-
           and trans-)
           \item list containing the reference proteome and the "fake" proteins
           added at the end with a patterned naming convention (sp|denovo_HF_fake_protein)
           made up of the concatenated potential hybrid peptides.}
}
\description{
This function takes in three mandatory inputs: (1) all denovo
candidates (2) database search results and (3) the corresponding proteome
fasta file. The function's role is to extract high confidence de novo peptides
and to search for their existence in the proteome, whether the entire peptide
sequence or its pair fragments (in one or two proteins).
}
\details{
This function is based on the published algorithm by Faridi et al.
(2018) for the identification and categorization of hybrid peptides. The function
described here adopts a slightly modified version of the algorithm for
computational efficiency. The function starts by extracting unassigned denovo
spectra where the Average Local Confidence (assigned by PEAKS software), is
equivalent to the ALC cutoff which is based on the median of the assigned spectra
(between denovo and database search). The sequences of all peptides are searched
against the reference proteome. If there is a hit then, then, the peptide sequence
within a spectrum group considered as being linear and each spectrum group is
is then filtered so as to keep the highest ALC-ranking spectra. Then, the rest
of the spectra (spectra that did not contain any sequence that had an entire
match in the proteome database) then undergo a "cutting" procedure where each
sequence yields n-2 sequences (with n being the length of the peptide. That is
if the peptide contains 9 amino acids i.e NTYASPRFK, then the sequence is cut
into a combination of 7 sequences of 2 fragment pairs each i.e fragment 1: NTY
and fragment 2: ASPRFK, etc).These are then searched in the proteome for hits
of both peptide fragments within a same protein, spectra in which sequences
have fragment pairs that match within a same protein, these are considerent
to be potentially cis-spliced. Potentially cis-spliced spectrum groups are then
filtered based on the highest ranking ALC. Spectrum groups not considered to be
potentially cis-spliced are further checked for potential trans-splicing. The
peptide sequences are cut again in the same fashion, however, this time peptide
fragment pairs are searched for matches in two proteins. Peptide sequences whose
fragment pairs match in 2 proteins are considerend to be potentially trans-spliced.
The same filtering for the highest ranking ALC within each peptide spectrum group.
The remaining spectra that were neither assigned as linear nor potentially spliced
(neither cis- nor trans-) are then discarded. The result is a list of spectra
along with their categorizations (Linear, potentially cis- and potentially trans-)
Potentially cis- and trans-spliced peptides are then concatenated and then broken into
several "fake" proteins and added to the bottom of the reference proteome. The
point of this last step is to create a merged proteome (consisting of the reference
proteome and the hybrid proteome) which  would be used for a second database
search. After the second database search the checknetmhcpan function or the
step2_wo_netMHCpan function can be used in order to obtain the final list
of potentially spliced peptides.
Article: Faridi P, Li C, Ramarathinam SH, Vivian JP, Illing PT, Mifsud NA,
Ayala R, Song J, Gearing LJ, Hertzog PJ, Ternette N, Rossjohn J, Croft NP,
Purcell AW. A subset of HLA-I peptides are not genomically templated: Evidence
for cis- and trans-spliced peptide ligands. Sci Immunol. 2018 Oct 12;3(28):eaar3947.
<doi: 10.1126/sciimmunol.aar3947>. PMID: 30315122.
}
\examples{
\dontrun{
 hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
 proteome, export = TRUE, output_dir)
 hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
 proteome)
 hybridFinderResult_list <- HybridFinder(denovo_candidates, db_search,
 proteome, export = FALSE)
}
}
\seealso{
\code{\link[seqinr]{read.fasta}},\code{\link[seqinr]{s2c}}
}
