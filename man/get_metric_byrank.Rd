\encoding{UTF-8}
\name{get_metric_byrank}
\alias{get_metric_byrank}
\title{Chemical metrics for taxa aggregated to a given rank}
\description{
Calculates a single chemical metric for taxa in each sample aggregated to a specified rank.
}

\usage{
  get_metric_byrank(RDP, map, refdb = "GTDB_220", taxon_AA = NULL,
    groups = NULL, zero_AA = NULL, metric = "Zc", rank = "genus")
}

\arguments{
  \item{RDP}{data frame, taxonomic abundances produced by \code{\link{read_RDP}} or \code{\link{ps_taxacounts}}}
  \item{map}{data frame, taxonomic mapping produced by \code{\link{map_taxa}}}
  \item{refdb}{character, name of reference database (\samp{GTDB_220} or \samp{RefSeq_206})}
  \item{taxon_AA}{data frame, amino acid compositions of taxa, used to bypass \code{refdb} specification}
  \item{groups}{list of indexing vectors, samples to be aggregated into groups}
  \item{zero_AA}{character, three-letter abbreviation(s) of amino acid(s) to assign zero counts for calculating chemical metrics}
  \item{metric}{character, chemical metric to calculate}
  \item{rank}{character, amino acid compositions of all lower-ranking taxa (to genus) are aggregated to this rank}
}

\details{

This function adds up amino acid compositions of taxa up to the specified rank and returns a data frame samples on the rows and taxa on the columns.
Because amino acid composition for genera have been precomputed from species-level genomes in a reference database, chemical metrics for genera are constant.
In contrast, chemical metrics for higher-level taxa is variable as they depend on the reference genomes as well as relative abundances of children taxa.

The value for \code{rank} should be one of \samp{rootrank}, \samp{domain}, \samp{phylum}, \samp{class}, \samp{order}, \samp{family}, or \samp{genus}.
For all ranks other than \samp{genus}, the amino acid compositions of all lower-ranking taxa are weighted by taxonomic abundance and summed in order to calculate the chemical metric at the specified rank.
If the rank is \samp{genus}, then no aggregation is done (because it is lowest-level rank available in the classifications), and the values of the \code{metric} for all genera in each sample are returned.
If the rank is \samp{rootrank}, then the results are equivalent to community reference proteomes (i.e., \code{\link{get_metrics}}).

The \code{RDP}, \code{map}, \samp{refdb}, and \code{groups} arguments are the same as described in \code{\link{get_metrics}}.
See \code{\link[canprot]{calc_metrics}} for available metrics.

}

\value{
A data frame of numeric values with row names corresponding to samples and column names corresponding to taxa.
}

\seealso{
\code{\link{get_metrics}}
}

\examples{
# Plot similar to Fig. 1 in Dick and Shock (2013)
# Read example dataset
RDPfile <- system.file("extdata/RDP-GTDB_220/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
# Get mapping to reference database
map <- map_taxa(RDP)
# Calculate phylum-level Zc
phylum_Zc <- get_metric_byrank(RDP, map, rank = "phylum")
# Keep phyla present in at least two samples
n_values <- colSums(!sapply(phylum_Zc, is.na))
phylum_Zc <- phylum_Zc[n_values > 2]
# Swap first two samples to get them in the right location
# (MG-RAST accession numbers for these samples are not in spatial order)
phylum_Zc <- phylum_Zc[c(2, 1, 3, 4, 5), ]
matplot(phylum_Zc, type = "b", xlab = "Sampling site (hot -> cool)", ylab = "Zc")
title("Phylum-level Zc at Bison Pool hot spring")
}

\references{

Dick JM, Shock E. 2013. A metastable equilibrium model for the relative abundances of microbial phyla in a hot spring. \emph{PLOS One} \bold{8}: e72395. \doi{10.1371/journal.pone.0072395}

}
