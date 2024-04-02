# Read and sum amino acid composition of all proteins for each genome
# 20221022 r207
# 20240318 r214
genome_AA <- function() {

  # Read the list of genome files for species representatives
  # wget "https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/gtdb_proteins_aa_reps_r214.tar.gz"
  bacfiles <- dir("214.1/genomic_files_reps/protein_faa_reps/bacteria/", full.names = TRUE)
  arcfiles <- dir("214.1/genomic_files_reps/protein_faa_reps/archaea/", full.names = TRUE)
  files <- c(bacfiles, arcfiles)

  # Loop over FASTA files (one for each genome)
  ifile <- seq_along(files)
  aa <- lapply(ifile, function(i) {
    # Print progress message
    if(i %% 100 == 0) print(i)
    # Read amino acid composition 
    aa <- suppressMessages(canprot::read_fasta(files[i]))
    # Sum amino acid composition
    canprot::sum_aa(aa)
  })
  aa <- do.call(rbind, aa)

  # Put in full genome names (with version suffix .1, .2, etc.)
  genome <- unlist(strsplit(basename(files[ifile]), "_protein.faa"))
  aa$organism <- genome

  # Put in species names 20240104
  taxonomy <- read.csv("taxonomy.csv")
  aa$ref <- taxonomy$species

  # Save result
  write.csv(aa, "genome_AA.csv", row.names = FALSE, quote = FALSE)

}
