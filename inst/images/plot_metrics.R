# Script for making plot_metrics.png (adapted from ?plot_metrics)
# 20230917 jmd

par(mar = c(4, 4, 3, 2))
par(cex = 1.5)

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
canprot::add_hull(pm$Zc[ilo], pm$nH2O[ilo],
  col = adjustcolor("green3", alpha.f = 0.3), border = NA)
canprot::add_hull(pm$Zc[ihi], pm$nH2O[ihi],
  col = adjustcolor("blue", alpha.f = 0.3), border = NA)

savePlot("plot_metrics.png")
