\encoding{UTF-8}
\name{get_metrics}
\alias{get_metrics}
\title{Calculate chemical metrics of community reference proteomes}
\description{
Combines RDP classifications with RefSeq reference proteomes to get amino acid compositions of community reference proteomes.
Amino acid compositions are used to calculate chemical metrics (\ZC, \nO2, and \nH2O).
}

\usage{
  get_metrics(RDP, map, refdb = "RefSeq", taxon_AA = NULL, groups = NULL,
    return_AA = FALSE, zero_AA = NULL)
}

\arguments{
  \item{RDP}{data frame, output from \code{\link{read_RDP}}}
  \item{map}{data frame, output from \code{\link{map_taxa}}}
  \item{refdb}{character, name of reference database (\samp{RefSeq} or \samp{GTDB})}
  \item{taxon_AA}{data frame, content of \code{extdata/RefSeq/taxon_AA.csv.xz} or \code{extdata/GTDB/taxon_AA.csv.xz} (used for speed-up)}
  \item{groups}{list, indexing vectors to identify groups of samples whose taxonomic classifications are aggregated before calculating amino acid composition and chemical metrics}
  \item{return_AA}{logical, return the amino acid composition for each sample instead of the chemical metrics?}
  \item{zero_AA}{character, set counts of the specified amino acids to zero for calculating chemical metrics (uses three-letter abbreviations).}
}

\details{

\code{get_metrics} calculates the carbon oxidation state (\ZC), stoichiometric oxidation state (\nO2), and stoichiometric hydration state (\nH2O) for the community reference proteome in each sample.
The community reference proteome is computed as the weighted mean of amino acid compositions of reference proteomes for taxa, where weights are the RDP classification counts.
Then, \ZC, \nO2, and \nH2O are calculated from the amino acid composition of the community reference proteome with the \code{\link[canprot]{ZCAA}}, \code{\link[canprot]{O2AA}}, and \code{\link[canprot]{H2OAA}} functions.

}

\seealso{
\code{\link{plot_metrics}} for plotting chemical metrics.
}

\examples{
## First two examples are for RDP Classifier with default training set
## and mapping to NCBI taxonomy with RefSeq reference proteomes

# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP)
# This is a data frame with 14 rows and Run, ZC, nO2, and nH2O columns
(metrics <- get_metrics(RDP, map))

# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
# Create list with metadata and metrics in same sample order
mdat <- get_metadata(mdatfile, metrics)
# Calculate metrics for aggregated samples of Archaea and Bacteria
groups <- list(A = mdat$metadata$domain == "Archaea", B = mdat$metadata$domain == "Bacteria")
# This is a data frame with 2 rows and group, ZC, nO2, and nH2O columns
get_metrics(RDP, map, groups = groups)

## Third example is for RDP Classifer with GTDB SSU training set
## and mapping to GTDB taxonomy with GTDB reference proteomes

RDPfile.GTDB <- system.file("extdata/RDP-GTDB/BGPF13.tab.xz", package = "chem16S")
RDP.GTDB <- read_RDP(RDPfile.GTDB)
map.GTDB <- map_taxa(RDP.GTDB)
metrics.GTDB <- get_metrics(RDP.GTDB, map.GTDB)

# Plot ZC from GTDB vs RefSeq
xylim <- range(metrics$ZC, metrics.GTDB$ZC)
plot(metrics$ZC, metrics.GTDB$ZC, xlim = xylim, ylim = xylim, type = "n")
lines(xylim, xylim, lty = 2, col = 8)
points(metrics$ZC, metrics.GTDB$ZC, pch = mdat$metadata$pch, bg = mdat$metadata$col)
md.leg <- mdat$metadata[1:2, ]
legend("bottomright", md.leg$domain, pch = md.leg$pch, pt.bg = md.leg$col)
title(quote(italic(Z)[C]~"from GTDB vs RefSeq"))

## Suppose we want to exclude tryptophan, tyrosine, and phenylalanine
## from the calculation of chemical metrics
get_metrics(RDP.GTDB, map.GTDB, zero_AA = c("Trp", "Tyr", "Phe"))
}
