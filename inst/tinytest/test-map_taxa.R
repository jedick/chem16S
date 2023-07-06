# 20220505

file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file, quiet = TRUE)
map <- map_taxa(RDP, refdb = "RefSeq", quiet = TRUE)

info <- "NA mappings are listed in unmapped_groups attribute"
expect_equal(sum(is.na(map)), length(attr(map, "unmapped_groups")), info = info)

info <- "Percentage of NA mapping is consistent with classification counts"
allcounts <- sum(RDP[, -(1:4)])
NAcounts <- sum(RDP[is.na(map), -(1:4)])
expect_equal(sum(attr(map, "unmapped_percent")), NAcounts / allcounts * 100, info = info)

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
expect_equal(sum(samenames, na.rm = TRUE), 53, info = info)
# Some different names are manually mapped
expect_equal(RDPname[sapply(samenames, isFALSE)], c("genus_GpIIa", "class_Cyanobacteria"), info = info)
expect_equal(NCBIname[sapply(samenames, isFALSE)], c("genus_Synechococcus", "phylum_Cyanobacteria"), info = info)

# 20230706

info <- "First example from ?map_taxa produces expected output"
file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file, quiet = TRUE)
expect_message(map <- map_taxa(RDP, refdb = "RefSeq"), "class_Cyanobacteria --> phylum_Cyanobacteria \\(12.3%\\)", info = info)
expect_equal(length(attributes(map)$unmapped_percent), 8, info = info)
expect_equal(round(sum(attributes(map)$unmapped_percent), 2), 24.42)

info <- "Sanity check that ranks are conserved, except for some groups"
# This dataset requires a fairly complex manual mapping
file <- system.file("extdata/RDP/HLA+16.tab.xz", package = "chem16S")
RDP <- read_RDP(file)
map <- map_taxa(RDP, refdb = "RefSeq")
# The following lines are adapted from get_metrics()
# Exclude NA mappings
RDP <- RDP[!is.na(map), ]
map <- na.omit(map)
# Get amino acid compositions of taxa compiled from RefSeq
refdb <- "RefSeq"
AApath <- file.path("extdata", refdb, "taxon_AA.csv.xz")
AAfile <- system.file(AApath, package = "chem16S")
taxon_AA <- read.csv(AAfile, as.is = TRUE)
# Apply the mapping
taxon_AA <- taxon_AA[map, ]
# Check that most of the ranks are unchanged in mapping
equalrank <- RDP$rank == taxon_AA$protein
expect_equal(sum(equalrank), 361, info = info)
# Check the names of the few mappings that cross ranks
icrossrank <- RDP$rank != taxon_AA$protein
expect_equal(RDP$name[icrossrank], c("Actinobacteria", "Spartobacteria_genera_incertae_sedis", "Subdivision3_genera_incertae_sedis", "Cyanobacteria"), info = info)
expect_equal(taxon_AA$organism[icrossrank], c("Actinobacteria", "Spartobacteria", "Verrucomicrobia subdivision 3", "Cyanobacteria"), info = info)
