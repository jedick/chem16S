# Get amino acid compositions for bacterial and archaeal taxa in GTDB
# 20221015 jmd

# Read and sum amino acid composition of all proteins for each genome 20221022
genome_AA <- function() {

  # The protein_faa_reps directory was found in this archive:
  # https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_proteins_aa_reps_r207.tar.gz
  bacfiles <- dir("207.0/genomic_files_reps/protein_faa_reps/bacteria/", full.names = TRUE)
  arcfiles <- dir("207.0/genomic_files_reps/protein_faa_reps/archaea/", full.names = TRUE)
  files <- c(bacfiles, arcfiles)

  # Loop over FASTA files (one for each genome)
  ifile <- seq_along(files)
  aa <- lapply(ifile, function(i) {
    # Print progress message
    if(i %% 100 == 0) print(i)
    # Read amino acid composition 
    aa <- suppressMessages(read.fasta(files[i]))
    # Sum amino acid composition
    aasum(aa)
  })
  aa <- do.call(rbind, aa)

  # Put in full genome names (with version suffix .1, .2, etc.)
  genome <- unlist(strsplit(basename(files[ifile]), "_protein.faa"))
  aa$organism <- genome

  # Save result
  write.csv(aa, "genome_AA.csv", row.names = FALSE, quote = FALSE)

}

# Combine amino acid composition of genomes for genus and higher levels 20221015
taxon_AA <- function() {

  # Read the summed amino acid compositions of all proteins for each genome
  genaa <- read.csv("genome_AA.csv")
  # Keep only genomes with at least 500 proteins
  genaa <- genaa[genaa$chains >= 500, ]
  # Normalize by number of proteins
  genaa[, 5:25] <- genaa[, 5:25] / genaa$chains

  # Read the GTDB taxonomy
  bactax <- read.table("207.0/bac120_taxonomy_r207.tsv.gz", sep = "\t")
  arctax <- read.table("207.0/ar53_taxonomy_r207.tsv.gz", sep = "\t")
  taxonomy <- rbind(bactax, arctax)

  # Match genomes to taxonomy
  itax <- match(genaa$organism, taxonomy[, 1])
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
    # Create blank data frame of amino acid composition (requires CHNOSZ > 1.4.3)
    aa <- seq2aa("", "")[rep(1, length(utaxa)), ]
    # Put in protein (rank) and organism (taxon) names
    aa$protein <- rank
    aa$organism <- utaxa
    # Loop over taxa
    for(i in 1:length(utaxa)) {
      # Sum amino acid compositions for all genomes in this taxon
      itaxa <- taxa == utaxa[i]
      aa[i, 5:25] <- aasum(genaa[itaxa, ])[, 5:25]
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
