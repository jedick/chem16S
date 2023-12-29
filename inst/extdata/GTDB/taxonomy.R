# Make taxonomy file for genomes 20231229
taxonomy <- function() {

  # Read the summed amino acid compositions of all proteins for each genome
  genaa <- read.csv("genome_AA.csv")
  # Read the GTDB taxonomy
  bactax <- read.table("207.0/bac120_taxonomy_r207.tsv.gz", sep = "\t")
  arctax <- read.table("207.0/ar53_taxonomy_r207.tsv.gz", sep = "\t")
  GTDBtax <- rbind(bactax, arctax)

  # Match genomes to taxonomy
  itax <- match(genaa$organism, GTDBtax[, 1])
  myGTDBtax <- GTDBtax[itax, ]
  # Get taxon names
  names <- strsplit(myGTDBtax[, 2], ";")
  # Remove d__, p__, c__, o__, f__, g__, s__ labels
  names <- lapply(names, function(x) gsub("^.__", "", x))
  taxonomy <- data.frame(
    genome = genaa$organism,
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
