# Combine amino acid composition of genomes for genus and higher levels 20221015
taxon_AA <- function() {

  # Read the summed amino acid compositions of all proteins for each genome
  genome_AA <- read.csv("genome_AA.csv")
  # Read the taxonomy 20231229
  taxonomy <- read.csv("taxonomy.csv")
  # Add a "rootrank" column 20241123
  taxonomy$rootrank = "Root"
  # Normalize by number of proteins
  genome_AA[, 5:25] <- genome_AA[, 5:25] / genome_AA$chains

  # Loop over taxonomic ranks
  ranks <- c("rootrank", "domain", "phylum", "class", "order", "family", "genus")
  aa <- lapply(ranks, function(rank) {
    # Names of all taxa at this rank
    taxa <- taxonomy[, rank]
    # Names of unique taxa
    utaxa <- unique(taxa)
    # Print rank and number of unique taxa
    print(paste(rank, length(utaxa)))
    # Create blank data frame of amino acid composition
    aa0 <- structure(list(protein = NA_character_, organism = NA_character_,
      ref = NA, abbrv = NA, chains = 0L, Ala = 0L, Cys = 0L, Asp = 0L,
      Glu = 0L, Phe = 0L, Gly = 0L, His = 0L, Ile = 0L, Lys = 0L,
      Leu = 0L, Met = 0L, Asn = 0L, Pro = 0L, Gln = 0L, Arg = 0L,
      Ser = 0L, Thr = 0L, Val = 0L, Trp = 0L, Tyr = 0L), row.names = "1", class = "data.frame")
    aa <- aa0[rep(1, length(utaxa)), ]
    # Put in protein (rank) and organism (taxon) names
    aa$protein <- rank
    aa$organism <- utaxa
    # Loop over taxa
    for(i in 1:length(utaxa)) {
      # Sum amino acid compositions for all genomes in this taxon
      itaxa <- taxa == utaxa[i]
      aa[i, 5:25] <- canprot::sum_aa(genome_AA[itaxa, ])[, 5:25]
    }
    # Normalize by number of genomes (put the number in 'ref' column)
    aa$ref <- aa$chains
    aa[, 5:25] <- aa[, 5:25] / aa$chains
    # For ranks below rootrank, put the higher-level rank in the 'abbrv' column
    irank <- match(rank, ranks)
    if(irank > 1) {
      uprank <- ranks[irank - 1]
      uptaxa <- taxonomy[, uprank]
      itaxa <- match(utaxa, taxa)
      aa$abbrv <- uptaxa[itaxa]
    }
    aa
  })

  aa <- do.call(rbind, aa)
  # Round the values
  aa[, 5:25] <- round(aa[, 5:25], 2)
  write.csv(aa, "taxon_AA.csv", row.names = FALSE, quote = FALSE)

}
