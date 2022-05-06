# 20220505

info <- "Calculated ZC increases along Bison Pool outflow channel"
# Get metrics for the Bison Pool dataset
RDPfile <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- mapRDP(RDP)
metrics <- getmetrics(RDP, map)
# Read the metadata file to put metrics in sample order
mdatfile <- system.file("extdata/metadata/SMS+12.csv", package = "chem16S")
mdat <- getmdat_example(mdatfile, metrics)
expect_true(all(diff(mdat$metrics$ZC) > 0))

# 20220506

info <- "Grouping samples puts names into output"
# Get chemical metrics for all samples in a dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- readRDP(RDPfile)
map <- mapRDP(RDP)
metrics <- getmetrics(RDP, map)
# Read the metadata file
mdatfile <- system.file("extdata/metadata/BGPF13.csv", package = "chem16S")
mdat <- getmdat_example(mdatfile)
# Make sure samples are in same order
stopifnot(mdat$Run == metrics$Run)
# Calculate metrics for aggregated samples of Archaea and Bacteria
groups0 <- list(mdat$cohort == "Archaea", mdat$cohort == "Bacteria")
met0 <- getmetrics(RDP, map, groups = groups0)
groups1 <- list(A = mdat$cohort == "Archaea", B = mdat$cohort == "Bacteria")
met1 <- getmetrics(RDP, map, groups = groups1)
expect_equal(met0$sample, c(1, 2))
expect_equal(met1$sample, c("A", "B"))
