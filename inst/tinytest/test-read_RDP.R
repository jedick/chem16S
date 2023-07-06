# 20220505

file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file, quiet = TRUE)

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
RDP.phylum <- read_RDP(file, lowest.level = "phylum", quiet = TRUE)
RDP.phylum.counts <- colSums(RDP.phylum[, -(1:4)])

# Test 1: column names are the same
expect_identical(colnames(RDP), colnames(RDP.phylum), info = info)
# Test 2: classification counts are the same
expect_identical(RDP.counts, RDP.phylum.counts, info = info)
# Test 3: only phyla are present in truncated result
expect_identical(unique(RDP.phylum$rank), "phylum", info = info)
