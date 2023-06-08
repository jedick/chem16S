\encoding{UTF-8}
\name{getmdat_example}
\alias{getmdat_example}
\title{Template function for metadata}
\description{
Simple function that processes metadata files.
}

\usage{
  getmdat_example(file, metrics)
}

\arguments{
  \item{file}{character, metadata file (CSV format)}
  \item{metrics}{data frame, output of \code{getmetrics}}
}

\details{
This function is provided as a simple example of reading metadata about 16S rRNA datasets.
The basic metadata for all datasets include \samp{name} (study name), \samp{Run} (SRA or other run ID), \samp{BioSample} (SRA BioSample), and \samp{Sample} (sample description).
Each data set contains environmental data, but the specific variables differ between datasets.
Therefore, this function processes the environmental data for different datasets and appends \samp{pch} and \samp{col} columns that can be used for visualization and identifying trends between groups of samples.

Samples are identified by the \code{Run} column (case-sensitive) in both \code{file} and \code{metrics}.
It is possible that metrics for one or more samples are not available (such as the number of RDP classifications is below \code{mincount} in \code{\link{getmetrics}}).
The function therefore removes metadata for unavailable samples.
Then, the samples in \code{metrics} are placed in the same order as they appear in metadata.
Ordering the samples is important because the columns in RDP Classifier output are not necessarily in the same order as the metadata.

In actual projects, customized getmdat_* functions can be used for more complex processing tasks.
}

\value{
  A list with two components: \code{metadata} and \code{metrics}.
}

\examples{
# Get metrics for the Bison Pool dataset (Swingley et al., 2012)
RDPfile <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- map_taxa(RDP)
metrics <- getmetrics(RDP, map)
# Read the metadata file
mdatfile <- system.file("extdata/metadata/SMS+12.csv", package = "chem16S")
mdat <- getmdat_example(mdatfile, metrics)

# Print sample names
mdat$metadata$Sample # 1 2 3 4 5
# Plot ZC of community reference proteomes at sites 1-5
Site <- mdat$metadata$Sample
ZC <- mdat$metrics$ZC
plot(Site, ZC, ylab = "Carbon oxidation state", type = "b")

# Find the columns of the RDP file for samples 1 to 5
match(mdat$metadata$Run, colnames(RDP)) # 6 5 7 8 9
# NOTE: for this dataset, MG-RAST IDs are in reversed order for samples 1 and 2
}
