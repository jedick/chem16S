# 20220505

info <- "Sample metrics are in same order as metadata"
# Get metrics for the Bison Pool dataset
RDPfile <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- map_taxa(RDP)
metrics <- getmetrics(RDP, map)
# Get metadata for the Bison Pool dataset
mdatfile <- system.file("extdata/metadata/SMS+12.csv", package = "chem16S")
mdat <- getmdat_example(mdatfile, metrics)
expect_true(all(mdat$metadata$Run == mdat$metrics$Run))
# This isn't a trivial test, because MG-RAST IDs are not in the sample order
expect_false(all(mdat$metadata$Run == metrics$Run))
