\encoding{UTF-8}
\name{read_RDP}
\alias{read_RDP}
\title{Read and filter RDP Classifier output}
\description{
Reads RDP Classifier output and optionally extracts lineages, removes classifications below a certain count, or truncates classifications to lowest taxonomic level.
}

\usage{
  read_RDP(file, lineage = NULL, mincount = 200,
    lowest.level = NULL, drop.groups = FALSE, quiet = FALSE)
}

\arguments{
  \item{file}{character, file name}
  \item{lineage}{character, regular expression for filtering lineages}
  \item{mincount}{integer, samples with less than this number of RDP classifications are excluded}
  \item{lowest.level}{character, truncate classifications at this taxonomic level}
  \item{drop.groups}{logical, remove some taxonomic groups deemed to be uninformative?}
  \item{quiet}{logical, suppress printed messages?}
}

\details{
\code{read_RDP} reads and filters RDP results for all samples in a study.
The input file should be created by RDP Classifer using the -h option to create a hierarchical listing.
Classifications for multiple samples can be combined into a single file using the merge-count command of RDP Classifier.

Only rows (lineages) with count > 0 for at least one sample are retained.

If \code{drop.groups} is TRUE:
Sequences with classifications at only the root or domain level (Root, Archaea, or Bacteria) are omitted because they provide poor taxonomic resolution.
Sequences classified to the class Chloroplast or genera Chlorophyta or Bacillariophyta are also omitted because they have little correspondence with the NCBI taxonomy.
These actions were hard-coded in earlier versions of chem16S (based on using the RDP Classifier taxonomy) but were made optional with adoption of the GTDB.

Then, only columns (samples) with classification count >= \code{mincount} are retained.
All remaining sequences (those classified to genus or higher levels) can be used for mapping to the NCBI taxonomy.

The lineage text of the RDP Classifier looks like \dQuote{Root;rootrank;Archaea;domain;Diapherotrites; phylum;Diapherotrites Incertae Sedis AR10;genus;}, so you can use \code{lineage = "Archaea"} to select the archaeal classifications or \code{lineage = "genus"} to select genus-level classifications.

Use the \code{lowest.level} argument to truncate the classifications to a level higher than genus.
This argument does not reduce the number of classifications, but only trims the RDP lineages to the specified level.
This may create duplicate lineages, for which the classification counts are summed, and only unique lineages are present in the returned data frame.

Change \code{quiet} to TRUE to suppress printing of messages about percentage classification to genus level, omitted sequences, and final range of total counts among all samples.
}

\value{
Data frame with columns inherited from \code{file} (from RDP output files: \code{taxid}, \code{rank}, \code{lineage}, \code{name}, and one column of classification counts for each sample).
If any sample has less than \code{mincount} total counts, that column is omitted from the result.
The rows in the result are subset from those in \code{file} by filtering operations.

NB: taxids in RDP files are not NCBI taxids.
}

\examples{
# An example file
file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")

# Default settings
r0 <- read_RDP(file)

# Suppress messages
r1 <- read_RDP(file, quiet = TRUE)

# Increase minimum count threshold
r2 <- read_RDP(file, mincount = 500)
# This should be 3 (i.e., 3 samples have less than 500 counts)
ncol(r1) - ncol(r2)

# Keep only Archaea
r3 <- read_RDP(file, lineage = "Archaea")
# This should be TRUE
all(grepl("Archaea", r3$lineage))

# Truncate taxonomy at phylum
r4 <- read_RDP(file, lowest.level = "phylum")
# This should be TRUE
all(sapply(strsplit(r4$lineage, ";"), tail, 1) == "phylum")
}
