# 20241119
# Read example dataset
RDPfile <- system.file("extdata/RDP/BGPF13.tab.xz", package = "chem16S")
RDP <- read_RDP(RDPfile)
map <- map_taxa(RDP, refdb = "RefSeq_206")

info <- "Get same results for get_metrics() and get_metric_byrank() at rootrank"
metric_byrank <- get_metric_byrank(RDP, map, refdb = "RefSeq_206", rank = "rootrank")
metrics <- get_metrics(RDP,  map, refdb = "RefSeq_206")
expect_equal(metric_byrank$Root, metrics$Zc, info = info)

info <- "Get same results using groups, a different metric, and zero_AA"
groups <- list(A = 1:7, B = 8:14)
metric_byrank <- get_metric_byrank(RDP, map, refdb = "RefSeq_206", groups = groups, zero_AA = "Met", metric = "nH2O", rank = "rootrank")
metrics <- get_metrics(RDP,  map, refdb = "RefSeq_206", groups = groups, zero_AA = "Met")
expect_equal(metric_byrank$Root, metrics$nH2O, info = info)
