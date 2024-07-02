\encoding{UTF-8}
\name{physeq}
\alias{physeq}
\alias{ps_taxacounts}
\alias{ps_metrics}
\alias{plot_ps_metrics}
\alias{plot_ps_metrics2}
\title{Process `phyloseq-class` objects to calculate chemical metrics}
\description{
Calculate chemical metrics of community reference proteomes from taxonomic abundances in \code{\link[phyloseq]{phyloseq-class}} objects and make plots of individual metrics or two metrics against each other.
}

\usage{
  ps_taxacounts(physeq, split = TRUE)
  ps_metrics(physeq, split = TRUE, refdb = "GTDB_214", quiet = FALSE, ...)
  plot_ps_metrics(physeq, metrics = c("Zc", "nO2", "nH2O"), x = "samples",
    color = NULL, shape = NULL, title = NULL,
    scales = "free_y", nrow = 1, sortby = NULL, ...)
  plot_ps_metrics2(physeq, metrics = c("Zc", "nH2O"), color = NULL,
    shape = NULL, title = NULL, refdb = "GTDB_214", quiet = FALSE)
}

\arguments{
  \item{physeq}{`phyloseq-class`, a \pkg{phyloseq} object containing an `otu_table` and `tax_table`}
  \item{split}{logical, calculate metrics for each sample (TRUE) or the entire set (FALSE)}
  \item{refdb}{character, name of reference database (passed to \code{\link{get_metrics}})}
  \item{quiet}{logical, suppress printed messages? (passed to \code{\link{map_taxa}})}
  \item{...}{additional arguments passed to \code{ps_metrics} and thence to \code{\link{get_metrics}}}
  \item{metrics}{character, chemical metrics to calculate}
  \item{x}{map this variable to the horizontal axis}
  \item{color}{map this variable to colors}
  \item{shape}{map this variable to shapes}
  \item{title}{main title for the plot}
  \item{scales}{\samp{free_y} or \samp{fixed} for different vertical scales or same scale in each panel}
  \item{nrow}{number of rows of panels (passed to \code{\link[ggplot2]{facet_wrap}})}
  \item{sortby}{subset of \code{metrics} argument, sort x-index by the mean values of these metrics}
}

\details{

\code{ps_taxacounts} returns a data frame with lowest-level classifications (from genus to phylum) and abundances for each OTU in the \code{physeq} object.

\code{ps_metrics} calculates chemical metrics of community reference proteomes from the lowest-level classifications after taxonomic mapping with \code{\link{map_taxa}}.

\code{plot_ps_metrics} plots individual chemical metrics or multiple chemical metrics, each in their own panel.
This function is modified from \code{\link[phyloseq]{plot_richness}}; see that help page for further details about graphics-related arguments.

\code{plot_ps_metrics2} plot two chemical metrics against each other.

The \code{metrics} argument is passed to \code{\link{get_metrics}} and selects the metrics to calculate.
\code{split}, \code{refdb}, and \code{quiet} can also be used as additional arguments in \code{plot_ps_metrics} and \code{plot_ps_metrics2}.

}

\value{
For \code{ps_taxacounts}, a data frame with samples in rows, and columns \code{taxid} (sample name taken from the OTU table obtained from \code{\link[phyloseq]{otu_table}(physeq)}), \code{lineage} (all taxonomic names from the taxonomy table obtained from \code{\link[phyloseq]{tax_table}(physeq)}, separated by semicolons), \code{name} (lowest-level taxonomic name), \code{rank} (rank of lowest-level taxon), followed by the abundances from the OTU table.

For \code{ps_metrics}, a data frame with samples in rows, and columns corresponding to each of the \code{metrics}.

For \code{plot_ps_metrics} and \code{plot_ps_metrics2}, ggplot objects.
}

\author{
The author of this Rd file and the \code{ps_taxacounts} and \code{ps_metrics} functions is Jeffrey Dick.
\code{plot_ps_metrics} is based on the \code{\link[phyloseq]{plot_richness}} function by Paul J. McMurdie and was modified from that function by Jeffrey Dick.
}

\examples{
# Example 1: Humboldt Sulfuretum (Fonseca et al., 2022)
file <- system.file("extdata/DADA2-GTDB_214/FEN+22/ps_FEN+22.rds", package = "chem16S")
physeq <- readRDS(file)
# Get lowest-level (to genus) classification for each OTU
taxacounts <- ps_taxacounts(physeq)
# Show numbers of assignments at each taxonomic level
table(taxacounts$rank)
# Map taxonomic names
map <- map_taxa(taxacounts, refdb = "GTDB_214")
# Taxonomy and reference proteomes are both from GTDB,
# so there are no unmapped classifications
sum(attr(map, "unmapped_percent")) # == 0
# Calculate chemical metrics
ps_metrics(physeq, refdb = "GTDB_214")

# Example 2: GlobalPatterns dataset from phyloseq
data(GlobalPatterns, package = "phyloseq")
# Plot metrics grouped by sample type and sorted by mean Zc;
# refdb = "RefSeq" uses manual mappings from the RDP to NCBI taxonomy
p <- plot_ps_metrics(GlobalPatterns, x = "SampleType",
                     sortby = "Zc", refdb = "RefSeq_206")
# Change orientation of x-axis labels
p + ggplot2::theme(axis.text.x = ggplot2::element_text(
  angle = 45, vjust = 1, hjust = 1))
}
