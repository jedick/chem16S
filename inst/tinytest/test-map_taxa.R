# 20220505

file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file)
map <- map_taxa(RDP, refdb = "RefSeq")

info <- "NA mappings are listed in unmapped_groups attribute"
expect_equal(sum(is.na(map)), length(attr(map, "unmapped_groups")))

info <- "Percentage of NA mapping is consistent with classification counts"
allcounts <- sum(RDP[, -(1:4)])
NAcounts <- sum(RDP[is.na(map), -(1:4)])
expect_equal(sum(attr(map, "unmapped_percent")), NAcounts / allcounts * 100)

info <- "For the test file, 53 RDP and NCBI names are identical and 2 are different"
# Read amino acid compositions of reference proteomes for genus and higher taxonomic groups from NCBI RefSeq`
AAfile <- system.file("extdata/RefSeq/taxon_AA.csv.xz", package = "chem16S")
AA <- read.csv(AAfile, as.is = TRUE)
# Make names by combining rank and name
NCBIname <- paste(AA$protein, AA$organism, sep = "_")[map]
RDPname <- paste(RDP$rank, RDP$name, sep = "_")
# Whether RDP and NCBI names are equal
samenames <- RDPname == NCBIname
# Identical names are automatically matched
expect_equal(sum(samenames, na.rm = TRUE), 53)
# Some different names are manually mapped
expect_equal(RDPname[sapply(samenames, isFALSE)], c("genus_GpIIa", "class_Cyanobacteria"))
expect_equal(NCBIname[sapply(samenames, isFALSE)], c("genus_Synechococcus", "phylum_Cyanobacteria"))
