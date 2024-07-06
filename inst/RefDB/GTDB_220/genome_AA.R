# Read and sum amino acid composition of all proteins for each genome
# 20221022 r207
# 20240318 r214
# 20240524 r220
genome_AA <- function() {

  # Read the list of genome files for species representatives
  # wget "https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/gtdb_proteins_aa_reps_r220.tar.gz"
  bacfiles <- dir("220.0/genomic_files_reps/protein_faa_reps/bacteria/", full.names = TRUE)
  arcfiles <- dir("220.0/genomic_files_reps/protein_faa_reps/archaea/", full.names = TRUE)
  files <- c(bacfiles, arcfiles)

  # Initialize data frame of amino acid composition
  ng <- length(files)
  aa <- data.frame(protein = character(ng), organism = character(ng),
    ref = character(ng), abbrv = character(ng), chains = numeric(ng),
    Ala = numeric(ng), Cys = numeric(ng), Asp = numeric(ng), Glu = numeric(ng),
    Phe = numeric(ng), Gly = numeric(ng), His = numeric(ng), Ile = numeric(ng),
    Lys = numeric(ng), Leu = numeric(ng), Met = numeric(ng), Asn = numeric(ng),
    Pro = numeric(ng), Gln = numeric(ng), Arg = numeric(ng), Ser = numeric(ng),
    Thr = numeric(ng), Val = numeric(ng), Trp = numeric(ng), Tyr = numeric(ng))

  # Loop over FASTA files (one for each genome)
  ifile <- seq_along(files)
  oldopt <- options(warn = 1)
  for(i in ifile) {
    # Print progress message
    if(i %% 100 == 0) print(i)
    # Read amino acid composition 
    myaa <- try(suppressMessages(canprot::read_fasta(files[i])), silent = TRUE)
    if(inherits(myaa, "try-error")) {
      print(paste0("Error reading file ", i, ": ", files[i]))
      next
    }
    # Sum amino acid composition
    aa[i, ] <- canprot::sum_aa(myaa)
  }
  options(oldopt)

  # Put in full genome names (with version suffix .1, .2, etc.)
  genome <- sapply(strsplit(basename(files[ifile]), "_protein.faa"), "[", 1)
  aa$organism <- genome

  # Put in species names 20240104
  taxonomy <- read.csv("taxonomy.csv")
  aa$ref <- taxonomy$species

  # Save result
  write.csv(aa, "genome_AA.csv", row.names = FALSE, quote = FALSE)

}
