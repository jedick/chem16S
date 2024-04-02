\encoding{UTF-8}
\name{plot_metrics}
\alias{plot_metrics}
\title{Plot chemical metrics of community reference proteomes}
\description{
Functions to plot chemical metrics of community reference proteomes.
}

\usage{
  plot_metrics(mdat,
    xvar = "Zc", yvar = "nH2O",
    xlim = NULL, ylim = NULL,
    xlab = chemlab(xvar), ylab = chemlab(yvar),
    plot.it = TRUE, add = FALSE, identify = FALSE,
    points = TRUE, lines = FALSE, title = TRUE,
    cex = 1, pt.open.col = 1, pch1 = 1, pch2 = 21,
    return = "data", extracolumn = NULL
  )
}

\arguments{
  \item{mdat}{list, output by \code{\link{get_metadata}}}
  \item{xvar}{character, name of x variable}
  \item{yvar}{character, name of y variable}
  \item{xlim}{numeric, x axis limits}
  \item{ylim}{numeric, y axis limits}
  \item{xlab}{x axis label}
  \item{ylab}{y axis label}
  \item{plot.it}{logical, make a plot?}
  \item{add}{logical, add to existing plot?}
  \item{identify}{logical, run \code{\link{identify}} for interactive identification of points?}
  \item{points}{logical, plot points?}
  \item{lines}{logical, plot lines?}
  \item{title}{character, plot title}
  \item{cex}{numeric, point size}
  \item{pt.open.col}{color of border for open point symbols (\code{pch > 20})}
  \item{pch1}{numeric, symbol for samples in group 1}
  \item{pch2}{numeric, symbol for samples in group 2}
  \item{return}{character, indicates whether to return \samp{data} values or group \samp{means}}
  \item{extracolumn}{character, the name of one or more extra columns (from \code{\link{get_metadata}}) to include in the output}
}

\details{

\code{plot_metrics} plots the values of \Zc and \nH2O (or other variables indicated by \code{xvar} and \code{yvar}) provided in \code{mdat$metrics}.
Symbol shape and color (\code{pch} and \code{col}) are taken from \code{mdat$metadata}.

If \code{pch1} and \code{pch2} are provided, then samples are classified into two groups according to the value of \code{mdat$metadata$pch}.
Mean values of the chemical metrics for each group are plotted with star-shaped symbols.

}

\value{
For \code{plot_metrics}, a data frame with columns for study name and Run IDs (\samp{name}, \samp{Run}), chemical metrics (taken from \code{mdat$metrics}), and symbols and colors for plotting points (\samp{pch}, \samp{col}).
}

\references{
Herlemann, D. P. R., Lundin, D., Andersson, A. F., Labrenz, M. and Jürgens, K. (2016) Phylogenetic signals of salinity and season in bacterial community composition across the salinity gradient of the Baltic Sea. \emph{Front. Microbiol.} \bold{7}, 1883. \doi{10.3389/fmicb.2016.01883}
}

\examples{
# Make a plot for the Baltic Sea salinity gradient
# (data from Herlemann et al., 2016)
RDPfile <- system.file("extdata/RDP/HLA+16.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq_206")
metrics <- get_metrics(RDP, map, refdb = "RefSeq_206")
mdatfile <- system.file("extdata/metadata/HLA+16.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)
pm <- plot_metrics(mdat)
# Add a legend
legend <- c("< 6 PSU", "6-20 PSU", "> 20 PSU")
pch <- c(24, 20, 21)
pt.bg <- c(3, NA, 4)
legend("bottomright", legend, pch = pch, col = 1, pt.bg = pt.bg, bg = "white")
# Classify samples with low and high salinity
ilo <- mdat$metadata$salinity < 6
ihi <- mdat$metadata$salinity > 20
# Add convex hulls
canprot::add_hull(pm$Zc[ilo], pm$nH2O[ilo],
  col = adjustcolor("green3", alpha.f = 0.3), border = NA)
canprot::add_hull(pm$Zc[ihi], pm$nH2O[ihi],
  col = adjustcolor("blue", alpha.f = 0.3), border = NA)

# Show points for all samples and larger star-shaped points
# for mean values of high- and low-salinity samples
plot_metrics(mdat, pch1 = 21, pch2 = 24)

# Plot nO2 instead of Zc
plot_metrics(mdat, xvar = "nO2")

# Make a plot for only Proteobacteria
RDP <- read_RDP(RDPfile, lineage = "Proteobacteria")
map <- map_taxa(RDP, refdb = "RefSeq_206")
metrics <- get_metrics(RDP, map, refdb = "RefSeq_206")
mdatfile <- system.file("extdata/metadata/HLA+16.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)
mdat$metadata$name <- paste(mdat$metadata$name, "(Proteobacteria)")
plot_metrics(mdat)
}
