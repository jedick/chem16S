\encoding{UTF-8}
\name{getmetrics}
\alias{getmetrics}
\title{Calculate chemical metrics from RDP classifications and RefSeq reference proteomes}
\description{
Combines RDP classifications with RefSeq reference proteomes to get amino acid compositions of community reference proteomes.
Amino acid compositions are used to calculate calculate chemical metrics (\ZC and \nH2O).
}

\usage{
  getmetrics(RDP, map, taxon_AA = NULL, groups = NULL)
}

\arguments{
  \item{RDP}{data frame, output from \code{\link{readRDP}}}
  \item{map}{data frame, output from \code{\link{mapRDP}}}
  \item{taxon_AA}{data frame, content of \code{extdata/RefSeq/taxon_AA.csv} (used for speed-up)}
  \item{groups}{list, indexing vectors to identify groups of samples whose taxonomic classifications are aggregated before calculating amino acid composition and chemical metrics}
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
# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- mapRDP(RDP)
metrics <- getmetrics(RDP, map)

# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
mdat <- getmdat_example(mdatfile)
# Make sure samples are in same order
stopifnot(mdat$Run == metrics$Run)

# Calculate metrics for aggregated samples of Archaea and Bacteria
groups <- list(A = mdat$cohort == "Archaea", B = mdat$cohort == "Bacteria")
getmetrics(RDP, map, groups = groups)
}
