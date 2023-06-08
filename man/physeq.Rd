\encoding{UTF-8}
\name{physeq}
\alias{physeq_taxacounts}
\alias{physeq_metrics}
\alias{plot_metrics}
\title{Calculate and plot chemical metrics from phyloseq objects}
\description{
Functions to process phyloseq objects for chemical metrics (\ZC, \nO2, and \nH2O) of community reference proteomes.
}

\usage{
  physeq_taxacounts(physeq, split = TRUE)
  physeq_metrics(physeq, split = TRUE, metrics = c("ZC", "nO2", "nH2O"), quiet = TRUE)
  plot_metrics(physeq, x = "samples", color = NULL, shape = NULL, title = NULL,
    scales = "free_y", nrow = 1, metrics = c("ZC", "nO2", "nH2O"), sortby = NULL)
}

\arguments{
  \item{physeq}{phyloseq-class, a phyloseq object containing an otu_table and tax_table}
  \item{split}{logical, calculate metrics for each sample (TRUE) or the entire set (FALSE)}
  \item{metrics}{character, names of metrics to calculate}
  \item{quiet}{logical, suppress printed messages? (passed to \code{\link{map_taxa}})}
  \item{x}{map this variable to horizontal axis}
  \item{color}{map this variable to different colors}
  \item{shape}{map this variable to different shapes}
  \item{title}{main title for the plot}
  \item{scales}{\samp{free_y} or \samp{fixed} for different vertical scales or same scale in each panel}
  \item{nrow}{number of rows of panels (passed to \code{\link{facet_wrap}})}
  \item{sortby}{subset of \code{metrics} argument, sort x-index by the mean values of these metrics}
}

\details{

\code{physeq_taxacounts} returns a data frame with lowest-level classifications (from genus to domain) and abundances for each OTU in the \code{physeq} object.
\code{physeq_metrics} calculates chemical metrics of community reference proteomes from the lowest-level classifications after mapping to the NCBI taxonomy with \code{\link{map_taxa}}.

\code{plot_metrics} creates a ggplot object for plotting the metrics.
This function is modified from \code{\link{plot_richness}}; see that help page for further details about graphics-related arguments.

}

\value{
For \code{physeq_taxacounts}, a data frame with samples in rows, and columns \code{taxid} (sample name taken from the OTU table obtained from \code{otu_table(physeq)}), \code{lineage} (all taxonomic names from the taxonomy table obtained from \code{tax_table(physeq)}, separated by semicolons), \code{name} (lowest-level taxonomic name), \code{rank} (rank of lowest-level taxon), followed by the abundances from the OTU table.

For \code{physeq_metrics}, a data frame with samples in rows, and columns corresponding to each of the \code{metrics}.

For \code{plot_metrics}, a ggplot object.
}

\author{
The author of this Rd file and the \code{physeq_taxacounts} and \code{physeq_metrics} functions is Jeffrey Dick.
The \code{\link{plot_richness}} function was written by Paul J. McMurdie and was modified by Jeffrey Dick to make \code{plot_metrics}.
}

\examples{
# Load the GlobalPatterns dataset from phyloseq
data(GlobalPatterns, package = "phyloseq")
# Get lowest-level (to genus) taxon for each OTU
taxacounts <- physeq_taxacounts(GlobalPatterns)
# Show numbers of assignments at each taxonomic level
table(taxacounts$rank)
# Calculate chemical metrics
# quiet = FALSE prints messages about taxonomic mapping
physeq_metrics(GlobalPatterns, quiet = FALSE)
# Plot metrics grouped by sample type
plot_metrics(GlobalPatterns, x = "SampleType")
# Plot metrics grouped by sample type and sorted by mean ZC
plot_metrics(GlobalPatterns, x = "SampleType", sortby = "ZC")
}
