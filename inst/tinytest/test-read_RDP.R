# 20220505

file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file, drop.groups = TRUE, quiet = TRUE)

info <- "Expected classifications are produced"
# These are five sites at Bison Pool from high to low temperature
# NOTE: mgm4946224.3 is site 1 (source pool) and mgm4946223.3 is site 2 (outflow channel)
# (i.e., sample location does not correspond to numerical order of first 2 IDs)
# Site 1, 93.3 °C
expect_equal(RDP$name[which.max(RDP$mgm4946224.3)], "Thermocrinis", info = info)
# Site 2, 79.4 °C
expect_equal(RDP$name[which.max(RDP$mgm4946223.3)], "Diapherotrites Incertae Sedis AR10", info = info)
# Site 3, 67.5 °C
expect_equal(RDP$name[which.max(RDP$mgm4946225.3)], "Thermus", info = info)
# Site 4, 65.3 °C
expect_equal(RDP$name[which.max(RDP$mgm4946226.3)], "Cyanobacteria", info = info)
# Site 5, 57.1 °C
expect_equal(RDP$name[which.max(RDP$mgm4946227.3)], "Enterobacter", info = info)

info <- "lowest.level truncates taxonomy but does not change classification counts"

# Default: classifications go down to genus level
RDP.counts <- colSums(RDP[, -(1:4)])
# Truncate classifications at phylum level
RDP.phylum <- read_RDP(file, lowest.level = "phylum", drop.groups = TRUE, quiet = TRUE)
RDP.phylum.counts <- colSums(RDP.phylum[, -(1:4)])

# Test 1: column names are the same
expect_identical(colnames(RDP), colnames(RDP.phylum), info = info)
# Test 2: classification counts are the same
expect_identical(RDP.counts, RDP.phylum.counts, info = info)
# Test 3: only phyla are present in truncated result
expect_identical(unique(RDP.phylum$rank), "phylum", info = info)

# 20241123
info <- "drop.groups = FALSE (the default) includes domain-level classifications"
RDP.phylum <- read_RDP(file, lowest.level = "phylum", quiet = TRUE)
expect_equal(unique(RDP.phylum$rank), c("phylum", "domain"), info = info)

# 20250321
info <- "'lineage' argument works as expected"
RDP.archaea <- read_RDP(file, lineage = "Archaea", quiet = TRUE)
expect_length(grep("Bacteria", RDP.archaea$lineage), 0, info = info)

info <- "Error if no samples have at least 'mincount' counts"
expect_error(read_RDP(file, lineage = "Euryarchaeota", quiet = TRUE), "No samples are left!", info = info)
info <- "More samples remain by lowering 'mincount'"
RDP.eury10 <- read_RDP(file, lineage = "Euryarchaeota", quiet = TRUE, mincount = 10)
RDP.eury1 <- read_RDP(file, lineage = "Euryarchaeota", quiet = TRUE, mincount = 1)
expect_true(ncol(RDP.eury1) > ncol(RDP.eury10), info = info)
