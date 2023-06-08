\encoding{UTF-8}
\name{readRDP}
\alias{readRDP}
\title{Read and filter RDP Classifier output}
\description{
Reads RDP Classifier output and optionally applies copy-number correction, extracts lineages, removes classifications below a certain count, or extracts lineages, or truncates classifications to lowest taxonomic level.
}

\usage{
  readRDP(file, lineage = NULL, mincount = 200,
    lowest.level = NULL, cn = FALSE, quiet = TRUE)
}

\arguments{
  \item{file}{character, file name}
  \item{lineage}{character, regular expression for filtering lineages}
  \item{mincount}{integer, samples with less than this number of RDP classifications are excluded}
  \item{lowest.level}{character, truncate classifications at this taxonomic level}
  \item{cn}{logical, adjust counts for 16S rRNA gene copy number?}
  \item{quiet}{logical, suppress printed messages?}
}

\details{
\code{readRDP} reads and filters RDP results for all samples in a study.
The input file should be created by RDP Classifer using the -h option to create a hierarchical listing.
Classifications for multiple samples can be combined into a single file using the merge-count command of RDP Classifier.

Only rows (lineages) with count > 0 for at least one sample are retained.
Sequences with classification results at only the root or domain level (Root, Archaea, or Bacteria) are omitted because they provide poor taxonomic resolution.
Sequences classified to the class Chloroplast or genera Chlorophyta or Bacillariophyta are also omitted because they have little correspondence with the NCBI taxonomy.
Then, only columns (samples) with classification count >= \code{mincount} are retained.
All remaining sequences (those classified to genus or higher levels) can be used for mapping to the NCBI taxonomy.

The lineage text of the RDP Classifier looks like \dQuote{Root;rootrank;Archaea;domain;Diapherotrites;phylum;Diapherotrites Incertae Sedis AR10;genus;}, so you can use \code{lineage = "Archaea"} to select the archaeal classifications or \code{lineage = "genus"} to select genus-level classifications.

Use the \code{lowest.level} argument to truncate the classifications to a level higher than genus.
This argument does not reduce the number of classifications, but only trims the RDP lineages to the specified level.
This may create duplicate lineages, for which the classification counts are summed, and only unique lineages are present in the returned data frame.

Change \code{quiet} to FALSE to enable printing of messages about percentage classification to genus level, omitted sequences, and final range of total counts among all samples.
}

\value{
  Data frame with columns inherited from \code{file} (for RDP output: taxid, rank, lineage, name, and one column of classification counts for each sample).
  If any sample has less than \code{mincount} counts, that column is omitted from the result.
  Output has rows from \code{file} that are selected or modified by filtering operations.

  NB: taxid is not the same as NCBI taxid.
}

\examples{
# An example file
file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")

# Default settings
r0 <- readRDP(file)

# Print a lot of messages
r1 <- readRDP(file, quiet = FALSE)

# Increase minimum count threshold
r2 <- readRDP(file, mincount = 500)
# This should be 3 (i.e., 3 samples have less than 500 counts)
ncol(r1) - ncol(r2)

# Apply copy-number correction
r3 <- readRDP(file, cn = TRUE)

# Keep only Archaea
r4 <- readRDP(file, lineage = "Archaea")
# This should be TRUE
all(grepl("Archaea", r4$lineage))

# Truncate taxonomy at phylum
r5 <- readRDP(file, lowest.level = "phylum")
# This should be TRUE
all(sapply(strsplit(r5$lineage, ";"), tail, 1) == "phylum")
}
