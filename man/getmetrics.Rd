\encoding{UTF-8}
\name{getmetrics}
\alias{getmetrics}
\title{Calculate chemical metrics from RDP classifications and RefSeq reference proteomes}
\description{
Combines RDP classifications with RefSeq reference proteomes to get amino acid compositions of community reference proteomes.
Amino acid compositions are used to calculate calculate chemical metrics (\ZC and \nH2O).
}

\usage{
  getmetrics(RDP, map, refdb = "RefSeq", taxon_AA = NULL, groups = NULL,
    return_AA = FALSE, zero_AA = NULL)
}

\arguments{
  \item{RDP}{data frame, output from \code{\link{readRDP}}}
  \item{map}{data frame, output from \code{\link{mapRDP}}}
  \item{refdb}{character, name of reference database (\samp{RefSeq} or \samp{GTDB})}
  \item{taxon_AA}{data frame, content of \code{extdata/RefSeq/taxon_AA.csv.xz} or \code{extdata/GTDB/taxon_AA.csv.xz} (used for speed-up)}
  \item{groups}{list, indexing vectors to identify groups of samples whose taxonomic classifications are aggregated before calculating amino acid composition and chemical metrics}
  \item{return_AA}{logical, return the amino acid composition for each sample instead of the chemical metrics?}
  \item{zero_AA}{character, set counts of the specified amino acids to zero for calculating chemical metrics (uses three-letter abbreviations).}
}

\details{

\code{getmetrics} calculates the carbon oxidation state (\ZC) and stoichiometric hydration state (\nH2O) for the community reference proteome in each sample.
The community reference proteome is computed as the weighted mean of amino acid compositions of reference proteomes for taxa, where weights are the RDP classification counts.
Then, \ZC and \nH2O are calculated from the amino acid composition of the community reference proteome with the \code{\link[canprot]{ZCAA}} and \code{\link[canprot]{H2OAA}} functions.

}

\seealso{
\code{\link{plotmet}} for plotting chemical metrics.
}

\examples{
## First two examples are for RDP Classifier with default training set
## and mapping to NCBI taxonomy with RefSeq reference proteomes

# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- mapRDP(RDP)
# This is a data frame with 14 rows and Run, nH2O, and ZC columns
(metrics <- getmetrics(RDP, map))

# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
# Create list with metadata and metrics in same sample order
mdat <- getmdat_example(mdatfile, metrics)
# Calculate metrics for aggregated samples of Archaea and Bacteria
groups <- list(A = mdat$metadata$domain == "Archaea", B = mdat$metadata$domain == "Bacteria")
# This is a data frame with 2 rows and group, nH2O, and ZC columns
getmetrics(RDP, map, groups = groups)

## Third example is for RDP Classifer with GTDB SSU training set
## and mapping to GTDB taxonomy with GTDB marker genes

RDPfile_GTDB <- system.file("extdata/RDP-GTDB/BGPF13.tab.xz", package = "chem16S")
RDP_GTDB <- readRDP(RDPfile_GTDB)
map_GTDB <- mapRDP(RDP_GTDB)
metrics_GTDB <- getmetrics(RDP_GTDB, map_GTDB)

# Plot ZC from GTDB vs RefSeq
xylim <- range(metrics$ZC, metrics_GTDB$ZC)
plot(metrics$ZC, metrics_GTDB$ZC, xlim = xylim, ylim = xylim, type = "n")
lines(xylim, xylim, lty = 2, col = 8)
points(metrics$ZC, metrics_GTDB$ZC, pch = mdat$metadata$pch, bg = mdat$metadata$col)
md.leg <- mdat$metadata[1:2, ]
legend("bottomright", md.leg$domain, pch = md.leg$pch, pt.bg = md.leg$col)
title(quote(italic(Z)[C]~"from GTDB vs RefSeq"))

## Suppose we want to exclude tryptophan, tyrosine, and phenylalanine
## from the calculation of chemical metrics
getmetrics(RDP_GTDB, map_GTDB, zero_AA = c("Trp", "Tyr", "Phe"))
}
