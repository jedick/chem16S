# Test added on 20250309

# Start NULL graphics device
pdf(NULL)
# Make a plot for the Baltic Sea salinity gradient
# (data from Herlemann et al., 2016)
RDPfile <- system.file("extdata/RDP/HLA+16.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq_206")
metrics <- get_metrics(RDP, map, refdb = "RefSeq_206")
mdatfile <- system.file("extdata/metadata/HLA+16.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)

info <- "plot_metrics() returns expected values for pch"
pm1 <- plot_metrics(mdat)
expect_equal(unique(pm1$pch), c(20, 24, 21), info = info)

info <- "plot_metrics() adds extra column of metadata"
pm2 <- plot_metrics(mdat, extracolumn = "type")
expect_equal(unique(pm2$type), c("Marine", "Mesohaline", "Oligohaline"), info = info)

info <- "plot_metrics() returns mean values of metrics"
pm3 <- plot_metrics(mdat, pch1 = 24, pch2 = 21, return = "means")
expect_equal(names(pm3), c("x1", "x2", "y1", "y2"), info = info)
# Mean nH2O is lower in marine than freshwater samples
expect_true(pm3$y2 < pm3$y1, info = info)

dev.off()
