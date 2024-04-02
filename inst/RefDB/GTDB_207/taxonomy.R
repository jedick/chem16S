# Make taxonomy file for genomes in GTDB 20231229
taxonomy <- function() {

  # The protein_faa_reps directory was found in this archive:
  # https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_proteins_aa_reps_r207.tar.gz
  bacfiles <- dir("207.0/genomic_files_reps/protein_faa_reps/bacteria/", full.names = TRUE)
  arcfiles <- dir("207.0/genomic_files_reps/protein_faa_reps/archaea/", full.names = TRUE)
  files <- c(bacfiles, arcfiles)
  # Get full genome names (with version suffix .1, .2, etc.)
  genome <- unlist(strsplit(basename(files), "_protein.faa"))

  # Read the GTDB taxonomy
  bactax <- read.table("207.0/bac120_taxonomy_r207.tsv.gz", sep = "\t")
  arctax <- read.table("207.0/ar53_taxonomy_r207.tsv.gz", sep = "\t")
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
