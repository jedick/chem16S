\encoding{UTF-8}
\name{get_metrics}
\alias{get_metrics}
\title{Calculate chemical metrics of community reference proteomes}
\description{
Combines RDP classifications with RefSeq reference proteomes to get amino acid compositions of community reference proteomes.
Amino acid compositions are used to calculate chemical metrics.
}

\usage{
  get_metrics(RDP, map, refdb = "GTDB", taxon_AA = NULL, groups = NULL,
    return_AA = FALSE, zero_AA = NULL, metrics = c("Zc", "nO2", "nH2O"))
}

\arguments{
  \item{RDP}{data frame, taxonomic abundances produced by \code{\link{read_RDP}} or \code{\link{ps_taxacounts}}}
  \item{map}{data frame, taxonomic mapping produced by \code{\link{map_taxa}}}
  \item{refdb}{character, name of reference database (\samp{GTDB} or \samp{RefSeq})}
  \item{taxon_AA}{data frame, contents of \file{taxon_AA.csv.xz} in \file{extdata/GTDB} or \file{extdata/RefSeq} (optionally used for speed-up)}
  \item{groups}{list of indexing vectors, samples to be aggregated into groups}
  \item{return_AA}{logical, return the amino acid composition for each sample instead of the chemical metrics?}
  \item{zero_AA}{character, three-letter abbreviation(s) of amino acid(s) to assign zero counts for calculating chemical metrics}
  \item{metrics}{character, the chemical metrics to calculate}
}

\details{

\code{get_metrics} calculates selected chemical metrics for the community reference proteome in each sample.
The community reference proteome is computed from the amino acid compositions of reference proteomes for taxa (obtained from the reference database in \code{refdb}), multiplied by taxonomic abundances given in \code{RDP}.
\code{RDP} may include results from the RDP Classifier (read using \code{\link{read_RDP}}) or derived from the OTU table of a \code{\link[phyloseq]{phyloseq-class}} object (see \code{\link{ps_taxacounts}}).
\code{map} defines the taxonomic mapping between \code{RDP} and \code{refdb}.
Then, chemical metrics are calculated from the amino acid composition of the community reference proteome.
The default chemical metrics are carbon oxidation state (\Zc), stoichiometric oxidation state (\nO2), and stoichiometric hydration state (\nH2O).
See \code{\link{calc_metrics}} for other available metrics.

\code{groups}, if given, is a list of one or more indexing vectors (with logical or numeric values) corresponding to samples whose taxonomic classifications are aggregated into groups before calculating amino acid compositions and chemical metrics.

}

\value{
A data frame with one row for each sample, corresponding to columns 5 and above of \code{RDP}.
The sample names are in the first column, which is named \code{Run} by default, or \code{group} if the \code{groups} argument is provided.
The remaining columns have numeric values and are named for each of the calculated \code{metrics}.
}

\seealso{
\code{\link{calc_metrics}}, \code{\link{plot_metrics}}
}

\examples{
## First two examples are for RDP Classifier with default training set
## and mapping to NCBI taxonomy with RefSeq reference proteomes

# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq")
# This is a data frame with 14 rows and Run, Zc, nO2, and nH2O columns
(metrics <- get_metrics(RDP, map, refdb = "RefSeq"))

# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
# Create list with metadata and metrics in same sample order
mdat <- get_metadata(mdatfile, metrics)
# Calculate metrics for aggregated samples of Archaea and Bacteria
groups <- list(A = mdat$metadata$domain == "Archaea",
  B = mdat$metadata$domain == "Bacteria")
# This is a data frame with 2 rows and group, Zc, nO2, and nH2O columns
get_metrics(RDP, map, refdb = "RefSeq", groups = groups)

## For the next example, classifications were made using the RDP Classifer
## with the GTDB SSU training set, and we use reference proteomes based on GTDB
## (this uses the default option of refdb = "GTDB").

RDPfile.GTDB <- system.file("extdata/RDP-GTDB/BGPF13.tab.xz", package = "chem16S")
RDP.GTDB <- read_RDP(RDPfile.GTDB)
map.GTDB <- map_taxa(RDP.GTDB)
metrics.GTDB <- get_metrics(RDP.GTDB, map.GTDB)

# Plot Zc from GTDB vs RefSeq
xylim <- range(metrics$Zc, metrics.GTDB$Zc)
plot(metrics$Zc, metrics.GTDB$Zc, xlim = xylim, ylim = xylim, type = "n")
lines(xylim, xylim, lty = 2, col = 8)
points(metrics$Zc, metrics.GTDB$Zc, pch = mdat$metadata$pch, bg = mdat$metadata$col)
md.leg <- mdat$metadata[1:2, ]
legend("bottomright", md.leg$domain, pch = md.leg$pch, pt.bg = md.leg$col)
title(quote(italic(Z)[C]~"from GTDB vs RefSeq"))

# To exclude tryptophan, tyrosine, and phenylalanine
# from the calculation of chemical metrics
get_metrics(RDP.GTDB, map.GTDB, zero_AA = c("Trp", "Tyr", "Phe"))
}
