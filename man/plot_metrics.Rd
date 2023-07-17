\encoding{UTF-8}
\name{plot_metrics}
\alias{plot_metrics}
\alias{add_hull}
\title{Plot chemical metrics of community reference proteomes}
\description{
Functions to plot chemical metrics of community reference proteomes.
}

\usage{
  plot_metrics(mdat, identify = FALSE, title = TRUE,
    xlim = NULL, ylim = NULL, plot.it = TRUE, points = TRUE, lines = FALSE,
    cex = 1, pch1 = 1, pch2 = 21, return = "data",
    extracolumn = NULL, add = FALSE, pt.open.col = 1,
    xlab = chemlab("Zc"), ylab = chemlab("nH2O"))
  add_hull(x, y, basecol, outline = FALSE, ...)
}

\arguments{
  \item{mdat}{list, output by \code{\link{get_metadata}}}
  \item{identify}{logical, run \code{\link{identify}} for interactive identification of points?}
  \item{title}{character, plot title}
  \item{xlim}{numeric, x axis limits}
  \item{ylim}{numeric, y axis limits}
  \item{plot.it}{logical, make a plot?}
  \item{points}{logical, plot points?}
  \item{lines}{logical, plot lines?}
  \item{cex}{numeric, point size}
  \item{pch1}{numeric, symbol for samples in group 1}
  \item{pch2}{numeric, symbol for samples in group 2}
  \item{return}{character, indicates whether to return \samp{data} values or group \samp{means}}
  \item{extracolumn}{character, the name of one or more extra columns (from \code{\link{get_metadata}}) to include in the output}
  \item{add}{logical, add to existing plot?}
  \item{pt.open.col}{color of border for open point symbols (\code{pch > 20})}
  \item{xlab}{x axis label}
  \item{ylab}{y axis label}
  \item{x}{x values}
  \item{y}{y values}
  \item{basecol}{color used for hull lines}
  \item{outline}{logical, draw hull lines?}
  \item{...}{arguments passed to \code{\link{polygon}}}
}

\details{

\code{plot_metrics} plots the \Zc and \nH2O values provided in \code{mdat$metrics}.
Symbol shape and color (\code{pch} and \code{col}) are taken from \code{mdat$metadata}.

If \code{pch1} and \code{pch2} are provided, then samples are classified into two groups according to the value of \code{mdat$metadata$pch}.
Mean values of the chemical metrics for each group are plotted with star-shaped symbols.

\code{add_hull} uses \code{\link{polygon}} to draw a convex hull around the points given in \code{x}, \code{y}.
The lines are drawn only if \code{outline} is TRUE; otherwise, the area of the hull is filled.
The fill color is a semi-transparent version of the color given in \code{basecol}.

}

\value{
For \code{plot_metrics}, a data frame with columns with study name and Run IDs (\samp{name}, \samp{Run}), plotting values (\samp{nH2O}, \samp{Zc}), and plotting symbols and colors (\samp{pch}, \samp{col}).
For \code{add_hull}, no return value (called for side effects).
}

\references{
Herlemann, D. P. R., Lundin, D., Andersson, A. F., Labrenz, M. and Jürgens, K. (2016) Phylogenetic signals of salinity and season in bacterial community composition across the salinity gradient of the Baltic Sea. \emph{Front. Microbiol.} \bold{7}, 1883. \doi{10.3389/fmicb.2016.01883}
}

\examples{
# Make a plot for the Baltic Sea salinity gradient
# (data from Herlemann et al., 2016)
RDPfile <- system.file("extdata/RDP/HLA+16.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq")
metrics <- get_metrics(RDP, map, refdb = "RefSeq")
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
add_hull(pm$Zc[ilo], pm$nH2O[ilo], "green3")
add_hull(pm$Zc[ihi], pm$nH2O[ihi], "blue")

# Show points for all samples and larger star-shaped points
# for mean values of high- and low-salinity samples
plot_metrics(mdat, pch1 = 21, pch2 = 24)

# Make a plot for only Proteobacteria
RDP <- read_RDP(RDPfile, lineage = "Proteobacteria")
map <- map_taxa(RDP, refdb = "RefSeq")
metrics <- get_metrics(RDP, map, refdb = "RefSeq")
mdatfile <- system.file("extdata/metadata/HLA+16.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)
mdat$metadata$name <- paste(mdat$metadata$name, "(Proteobacteria)")
plot_metrics(mdat)
}
