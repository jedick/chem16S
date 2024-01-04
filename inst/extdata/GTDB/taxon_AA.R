# Combine amino acid composition of genomes for genus and higher levels 20221015
taxon_AA <- function() {

  # Read the summed amino acid compositions of all proteins for each genome
  genome_AA <- read.csv("genome_AA.csv.xz")
  # Read the taxonomy 20231229
  taxonomy <- read.csv("taxonomy.csv")
  # Keep only genomes with at least 500 proteins
  i500 <- genome_AA$chains >= 500
  genome_AA <- genome_AA[i500, ]
  taxonomy <- taxonomy[i500, ]
  # Normalize by number of proteins
  genome_AA[, 5:25] <- genome_AA[, 5:25] / genome_AA$chains

  # Loop over taxonomic ranks
  ranks <- c("domain", "phylum", "class", "order", "family", "genus")
  aa <- lapply(ranks, function(rank) {
    # Names of all taxa at this rank
    taxa <- taxonomy[, rank]
    # Names of unique taxa
    utaxa <- unique(taxa)
    # Print rank and number of unique taxa
    print(paste(rank, length(utaxa)))
    # Create blank data frame of amino acid composition (requires CHNOSZ > 1.4.3)
    aa <- seq2aa("", "")[rep(1, length(utaxa)), ]
    # Put in protein (rank) and organism (taxon) names
    aa$protein <- rank
    aa$organism <- utaxa
    # Loop over taxa
    for(i in 1:length(utaxa)) {
      # Sum amino acid compositions for all genomes in this taxon
      itaxa <- taxa == utaxa[i]
      aa[i, 5:25] <- aasum(genome_AA[itaxa, ])[, 5:25]
    }
    # Normalize by number of genomes (put the number in 'ref' column)
    aa$ref <- aa$chains
    aa[, 5:25] <- aa[, 5:25] / aa$chains
    # For ranks below domain, put the higher-level rank in the 'abbrv' column
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
