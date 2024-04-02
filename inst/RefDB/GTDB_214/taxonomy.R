# Make taxonomy file for genomes in GTDB
# 20231229 r207
# 20240318 r214
taxonomy <- function() {

  # Read the list of genome files for species representatives
  # wget "https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/gtdb_proteins_aa_reps_r214.tar.gz"
  bacfiles <- dir("214.1/genomic_files_reps/protein_faa_reps/bacteria/", full.names = TRUE)
  arcfiles <- dir("214.1/genomic_files_reps/protein_faa_reps/archaea/", full.names = TRUE)
  files <- c(bacfiles, arcfiles)
  # Get full genome names (with version suffix .1, .2, etc.)
  genome <- sapply(strsplit(basename(files), "_protein.faa"), "[", 1)

  # Read the GTDB taxonomy
  # wget "https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv.gz"
  # wget "https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_taxonomy_r214.tsv.gz"
  bactax <- read.table("214.1/bac120_taxonomy_r214.tsv.gz", sep = "\t")
  arctax <- read.table("214.1/ar53_taxonomy_r214.tsv.gz", sep = "\t")
  GTDBtax <- rbind(bactax, arctax)

  # Match genomes to taxonomy
  itax <- match(genome, GTDBtax[, 1])
  myGTDBtax <- GTDBtax[itax, ]
  # Get taxon names
  names <- strsplit(myGTDBtax[, 2], ";")
  # Remove d__, p__, c__, o__, f__, g__, s__ labels
  names <- lapply(names, function(x) gsub("^.__", "", x))
  taxonomy <- data.frame(
    genome = genome,
    domain = sapply(names, "[", 1),
    phylum = sapply(names, "[", 2),
    class = sapply(names, "[", 3),
    order = sapply(names, "[", 4),
    family = sapply(names, "[", 5),
    genus = sapply(names, "[", 6),
    species = sapply(names, "[", 7)
  )
  write.csv(taxonomy, "taxonomy.csv", row.names = FALSE, quote = FALSE)

}
