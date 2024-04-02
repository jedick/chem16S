\encoding{UTF-8}
\name{get_metadata}
\alias{get_metadata}
\title{Template function for metadata}
\description{
Simple function that processes metadata files.
}

\usage{
  get_metadata(file, metrics)
}

\arguments{
  \item{file}{character, metadata file (CSV format)}
  \item{metrics}{data frame, output of \code{get_metrics}}
}

\details{
This function is provided as a simple example of reading metadata about 16S rRNA datasets.
The basic metadata for all datasets include \samp{name} (study name), \samp{Run} (SRA or other run ID), \samp{BioSample} (SRA BioSample), and \samp{Sample} (sample description).
Each data set contains environmental data, but the specific variables differ between datasets.
Therefore, this function processes the environmental data for different datasets and appends \samp{pch} and \samp{col} columns that can be used for visualization and identifying trends between groups of samples.

Samples are identified by the \code{Run} column (case-sensitive) in both \code{file} and \code{metrics}.
It is possible that metrics for one or more samples are not available (for example, if the number of RDP classifications is below \code{mincount} in \code{\link{get_metrics}}).
The function therefore removes metadata for unavailable samples.
Then, the samples in \code{metrics} are placed in the same order as they appear in metadata.
Ordering the samples is important because the columns in RDP Classifier output are not necessarily in the same order as the metadata.
}

\value{
If \code{metrics} is NULL, a data frame read from \code{file} with columns \samp{pch} and \samp{col} appended.
Otherwise, a list with two components, \code{metadata} and \code{metrics}.
In this list, \code{metrics} is the value from the \code{metrics} argument and \code{metadata} is a data frame read from \code{file} with columns \samp{pch} and \samp{col} appended and excluding any rows that correspond to unavailable samples in \code{metrics}.
}

\examples{
# Get metrics for the Bison Pool dataset (Swingley et al., 2012)
RDPfile <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq_206")
metrics <- get_metrics(RDP, map, refdb = "RefSeq_206")
# Read the metadata file
mdatfile <- system.file("extdata/metadata/SMS+12.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)

# Print sample names
mdat$metadata$Sample # 1 2 3 4 5
# Plot Zc of community reference proteomes at sites 1-5
Site <- mdat$metadata$Sample
Zc <- mdat$metrics$Zc
plot(Site, Zc, ylab = "Carbon oxidation state", type = "b")

# Find the columns of the RDP file for samples 1 to 5
match(mdat$metadata$Run, colnames(RDP)) # 6 5 7 8 9
# NOTE: for this dataset, samples 1 and 2 are in columns 6 and 5, respectively
}
