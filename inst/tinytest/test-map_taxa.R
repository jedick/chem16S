# 20220505

file <- system.file("extdata/RDP/SMS+12.tab.xz", package = "chem16S")
RDP <- read_RDP(file, quiet = TRUE)
map <- map_taxa(RDP, refdb = "RefSeq_206", quiet = TRUE)

info <- "NA mappings are listed in unmapped_groups attribute"
expect_equal(sum(is.na(map)), length(attr(map, "unmapped_groups")), info = info)

info <- "Percentage of NA mapping is consistent with classification counts"
allcounts <- sum(RDP[, -(1:4)])
NAcounts <- sum(RDP[is.na(map), -(1:4)])
expect_equal(sum(attr(map, "unmapped_percent")), NAcounts / allcounts * 100, info = info)

info <- "For the test file, 53 RDP and NCBI names are identical and 2 are different"
# Read amino acid compositions of reference proteomes for genus and higher taxonomic groups from NCBI RefSeq`
AAfile <- system.file("RefDB/RefSeq_206/taxon_AA.csv.xz", package = "chem16S")
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
RDP <- read_RDP(file, drop.groups = TRUE, quiet = TRUE)
expect_message(map <- map_taxa(RDP, refdb = "RefSeq_206"), "class_Cyanobacteria --> phylum_Cyanobacteria \\(12.3%\\)", info = info)
expect_equal(length(attributes(map)$unmapped_percent), 8, info = info)
expect_equal(round(sum(attributes(map)$unmapped_percent), 2), 24.42)

info <- "Sanity check that ranks are conserved, except for some groups"
# This dataset requires a fairly complex manual mapping
file <- system.file("extdata/RDP/HLA+16.tab.xz", package = "chem16S")
RDP <- read_RDP(file)
map <- map_taxa(RDP, refdb = "RefSeq_206")
# The following lines are adapted from get_metrics()
# Exclude NA mappings
RDP <- RDP[!is.na(map), ]
map <- na.omit(map)
# Get amino acid compositions of taxa compiled from RefSeq
refdb <- "RefSeq_206"
AApath <- file.path("RefDB", refdb, "taxon_AA.csv.xz")
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

# 20241130

## From https://data.gtdb.ecogenomic.org/releases/release220/220.0/RELEASE_NOTES.txt
## - Post-curation cycle, we identified updated spelling for 15 taxon names:
##    p__Calescibacterota (updated name: Calescibacteriota)
##    c__Brachyspirae (updated name: Brachyspiria)
##    c__Leptospirae (updated name: Leptospiria)
##    o__Ammonifexales (updated name: Ammonificales)
##    o__Exiguobacterales (updated name: Exiguobacteriales)
##    o__Hydrogenedentiales (updated name: Hydrogenedentales)
##    o__Phormidesmiales (updated name: Phormidesmidales)
##    f__Arcanobacteraceae (updated name: Arcanibacteraceae)
##    f__Acetonemaceae (updated name: Acetonemataceae)
##    f__Ethanoligenenaceae (updated name: Ethanoligenentaceae)
##    f__Exiguobacteraceae (updated name: Exiguobacteriaceae)
##    f__Geitlerinemaceae (updated name: Geitlerinemataceae)
##    f__Koribacteraceae (updated name: Korobacteraceae)
##    f__Phormidesmiaceae (updated name: Phormidesmidaceae)
##    f__Porisulfidaceae (updated name: Poriferisulfidaceae)
##   Note that the LPSN linkouts point to the correct updated names. We encourage users to use the updated names as these will appear in the next release.
## - Post-curation cycle, we discovered that two provisionally named families, Nitrincolaceae and Denitrovibrionaceae  have been validly named under the ICNP as Balneatricaceae and Geovibrionaceae, respectively.
##   We encourage users to use the validly published names as these will appear in the next release.

info <- "GTDB release 220 maps post-curation updated names"
# Get amino acid compositions of taxa compiled from GTDB
GTDB_220 = read.csv(system.file("RefDB/GTDB_220/taxon_AA.csv.xz", package = "chem16S"))
# Read phyloseq file generated with DADA2 pipeline and GTDB 16S rRNA sequences
psfile <- system.file("extdata/DADA2-GTDB_220/ZFZ+23/ps_ZFZ+23.rds", package = "chem16S")
ps <- readRDS(psfile)
taxa <- na.omit(unlist(data.frame(phyloseq::tax_table(ps))))
# We see that some deprecated names are present
oldnames <- c("Leptospirae", "Ammonifexales", "Phormidesmiales", "Exiguobacterales", 
"Hydrogenedentiales", "Phormidesmiaceae", "Exiguobacteraceae", 
"Koribacteraceae", "Acetonemaceae", "Nitrincolaceae")
expect_equal(setdiff(taxa, GTDB_220$organism), oldnames, info = info)
# map_taxa() takes care of mapping the deprecated names to current ones
taxacounts <- ps_taxacounts(ps)
expect_stdout(map <- map_taxa(taxacounts), "using these post-curation mapping\\(s\\)", info = info)
expect_false(any(is.na(map)), info = info)
