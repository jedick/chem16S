# Get amino acid compositions of bacterial and archaeal marker genes in GTDB
# 20221015 jmd

# Read the full GTDB taxonomy
bactax <- read.table("207.0/bac120_taxonomy_r207.tsv.gz", sep = "\t")
arctax <- read.table("207.0/ar53_taxonomy_r207.tsv.gz", sep = "\t")
taxonomy <- rbind(bactax, arctax)

# Get amino acid composition of concatenated marker genes for each genome 20221015
concat_marker <- function() {

  # Create an empty data frame of amino acid composition
  # Requires CHNOSZ > 1.4.3
  aa <- seq2aa("", "")[rep(1, nrow(taxonomy)), ]
  # Put in protein names
  nbac <- nrow(bactax)
  narc <- nrow(arctax)
  aa$protein[1:nbac] <- "bac120"
  aa$protein[(nbac + 1) : (nbac + narc)] <- "ar53"
  # Put in genome names
  aa$organism[1:nbac] <- bactax[, 1]
  aa$organism[(nbac + 1) : (nbac + narc)] <- arctax[, 1]

  # Loop over bacteria and archaea
  # These directories are unzipped from bac120_marker_genes_reps_r207.tar.gz and ar53_marker_genes_reps_r207.tar.gz
  for(bacarc in c("bac120_marker_genes_reps_r207", "ar53_marker_genes_reps_r207")) {

    # Loop over FASTA files (one for each marker gene)
    files <- dir(file.path("207.0", bacarc, "faa"), full.names = TRUE)
    for(file in files) {
      print(file)
      # Get amino acid compositions
      faa <- suppressMessages(read.fasta(file))
      # Match taxa to sequences
      ifaa <- match(taxonomy[, 1], faa$protein)
      faa <- faa[ifaa, ]
      # Zero out NAs
      faa[is.na(ifaa), 5:25] <- 0
      # Add amino acid composition to output
      aa[, 5:25] <- aa[, 5:25] + faa[, 5:25]
    }

  }

  # Remove genomes with 0 marker sequences
  aa <- aa[aa$chains != 0, ]
  # Save result
  write.csv(aa, "concat_marker_aa.csv", row.names = FALSE, quote = FALSE)

}

# Combine amino acid composition of genomes for genus and higher levels 20221015
taxon_AA <- function() {

  # Read the amino acid compositions of concatenated marker genes for each genome
  cmaa <- read.csv("concat_marker_aa.csv")
  # Keep only genomes with at least 80/120 (66%) markers for bacteria and 35/53 (66%) markers for archaea
  ibac <- cmaa$protein == "bac120" & cmaa$chains >= 80
  iarc <- cmaa$protein == "ar53" & cmaa$chains >= 35
  cmaa <- cmaa[ibac | iarc, ]
  # Normalize by number of marker genes
  cmaa[, 5:25] <- cmaa[, 5:25] / cmaa$chains

  # Match genomes with marker sequences to taxonomy
  itax <- match(cmaa$organism, taxonomy[, 1])
  mytaxonomy <- taxonomy[itax, ]
  # Get taxon names
  names <- strsplit(mytaxonomy[, 2], ";")
  # Remove d__, p__, c__, o__, f__, g__, s__ labels
  names <- lapply(names, function(x) gsub("^.__", "", x))
  domain <- sapply(names, "[", 1)
  phylum <- sapply(names, "[", 2)
  class <- sapply(names, "[", 3)
  order <- sapply(names, "[", 4)
  family <- sapply(names, "[", 5)
  genus <- sapply(names, "[", 6)
  species <- sapply(names, "[", 7)

  # Loop over taxonomic ranks
  ranks <- c("domain", "phylum", "class", "order", "family", "genus")
  aa <- lapply(ranks, function(rank) {
    # Names of all taxa at this rank
    taxa <- get(rank)
    # Names of unique taxa
    utaxa <- unique(taxa)
    # Print rank and number of unique taxa
    print(paste(rank, length(utaxa)))
    # Create blank data frame of amino acid composition
    aa <- seq2aa("", "")[rep(1, length(utaxa)), ]
    # Put in protein (rank) and organism (taxon) names
    aa$protein <- rank
    aa$organism <- utaxa
    # Loop over taxa
    for(i in 1:length(utaxa)) {
      # Sum amino acid compositions for marker genes of all genomes in this taxon
      itaxa <- taxa == utaxa[i]
      aa[i, 5:25] <- aasum(cmaa[itaxa, ])[, 5:25]
    }
    # Normalize by number of genomes (put the number in 'ref' column)
    aa$ref <- aa$chains
    aa[, 5:25] <- aa[, 5:25] / aa$chains
    # For ranks below domain, put the higher-level rank in the 'abbrv' column
    irank <- match(rank, ranks)
    if(irank > 1) {
      uprank <- ranks[irank - 1]
      uptaxa <- get(uprank)
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
