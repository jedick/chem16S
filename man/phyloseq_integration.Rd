\encoding{UTF-8}
\name{physeq}
\alias{ps_taxacounts}
\alias{ps_metrics}
\alias{plot_ps_metrics}
\alias{plot_ps_metrics2}
\title{Calculate and plot chemical metrics from phyloseq objects}
\description{
Functions to process phyloseq objects for chemical metrics (\Zc, \nO2, and \nH2O) of community reference proteomes.
}

\usage{
  ps_taxacounts(physeq, split = TRUE)
  ps_metrics(physeq, split = TRUE, metrics = c("Zc", "nO2", "nH2O"), quiet = TRUE, ...)
  plot_ps_metrics(physeq, x = "samples", color = NULL, shape = NULL, title = NULL,
    scales = "free_y", nrow = 1, metrics = c("Zc", "nO2", "nH2O"), sortby = NULL)
  plot_ps_metrics2(physeq, x = "Zc", y = "nH2O", color = NULL, shape = NULL, title = NULL)
}

\arguments{
  \item{physeq}{phyloseq-class, a phyloseq object containing an otu_table and tax_table}
  \item{split}{logical, calculate metrics for each sample (TRUE) or the entire set (FALSE)}
  \item{metrics}{character, names of metrics to calculate}
  \item{quiet}{logical, suppress printed messages? (passed to \code{\link{map_taxa}})}
  \item{...}{additional arguments for \code{\link{get_metrics}}}
  \item{x}{map this variable to the horizontal axis}
  \item{color}{map this variable to different colors}
  \item{shape}{map this variable to different shapes}
  \item{title}{main title for the plot}
  \item{scales}{\samp{free_y} or \samp{fixed} for different vertical scales or same scale in each panel}
  \item{nrow}{number of rows of panels (passed to \code{\link{facet_wrap}})}
  \item{sortby}{subset of \code{metrics} argument, sort x-index by the mean values of these metrics}
  \item{y}{map this variable to the vertical axis}
}

\details{

\code{ps_taxacounts} returns a data frame with lowest-level classifications (from genus to domain) and abundances for each OTU in the \code{physeq} object.
\code{ps_metrics} calculates chemical metrics of community reference proteomes from the lowest-level classifications after mapping to the NCBI taxonomy with \code{\link{map_taxa}}.

\code{plot_ps_metrics} plots individual chemical metrics or multiple chemical metrics, each in their own panel.
This function is modified from \code{\link{plot_richness}}; see that help page for further details about graphics-related arguments.

\code{plot_ps_metrics2} makes a cross-plot of two chemical metrics.

}

\value{
For \code{ps_taxacounts}, a data frame with samples in rows, and columns \code{taxid} (sample name taken from the OTU table obtained from \code{\link{otu_table}(physeq)}), \code{lineage} (all taxonomic names from the taxonomy table obtained from \code{\link{tax_table}(physeq)}, separated by semicolons), \code{name} (lowest-level taxonomic name), \code{rank} (rank of lowest-level taxon), followed by the abundances from the OTU table.

For \code{ps_metrics}, a data frame with samples in rows, and columns corresponding to each of the \code{metrics}.

For \code{plot_ps_metrics} and \code{plot_ps_metrics2}, ggplot objects.
}

\author{
The author of this Rd file and the \code{ps_taxacounts} and \code{ps_metrics} functions is Jeffrey Dick.
\code{plot_ps_metrics} is based on the \code{\link{plot_richness}} function by Paul J. McMurdie and was modified from that function by Jeffrey Dick.
}

\examples{
# Load the GlobalPatterns dataset from phyloseq
data(GlobalPatterns, package = "phyloseq")

\dontrun{
# Get lowest-level (to genus) taxon for each OTU
taxacounts <- ps_taxacounts(GlobalPatterns)
# Show numbers of assignments at each taxonomic level
table(taxacounts$rank)
# Map taxonomic names to NCBI taxonomy
map <- map_taxa(taxacounts)
# Show most prevalent unmapped groups
head(sort(table(attr(map, "unmapped_groups")), TRUE), 5)
# Calculate chemical metrics
# quiet = FALSE prints messages about taxonomic mapping
ps_metrics(GlobalPatterns, quiet = FALSE)
}

# Plot metrics grouped by sample type and sorted by mean Zc
plot_ps_metrics(GlobalPatterns, x = "SampleType", sortby = "Zc")
}
