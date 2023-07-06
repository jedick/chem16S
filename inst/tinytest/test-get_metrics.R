# 20220505
info <- "Calculated Zc increases along Bison Pool outflow channel"
# Get metrics for the Bison Pool dataset
RDPfile <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile, quiet = TRUE)
map <- map_taxa(RDP, refdb = "RefSeq", quiet = TRUE)
metrics <- get_metrics(RDP, map, refdb = "RefSeq")
# Read the metadata file to put metrics in sample order
mdatfile <- system.file("extdata/metadata/SMS+12.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)
expect_true(all(diff(mdat$metrics$Zc) > 0), info = info)

# 20220506
info <- "Grouping samples puts names into output"
# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile, quiet = TRUE)
map <- map_taxa(RDP, refdb = "RefSeq", quiet = TRUE)
metrics <- get_metrics(RDP, map, refdb = "RefSeq")
# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
mdat <- get_metadata(mdatfile, metrics)
# Calculate metrics for aggregated samples of Archaea and Bacteria
groups <- list(A = mdat$metadata$domain == "Archaea", B = mdat$metadata$domain == "Bacteria")
g.metrics <- get_metrics(RDP, map, groups = groups, refdb = "RefSeq")
expect_equal(g.metrics$group, c("A", "B"), info = info)

# 20221016
info <- "Bacteria are more oxidized than Archaea (ensures we don't have NA/NaN results)"
expect_true(diff(g.metrics$Zc) > 0, info = info)
