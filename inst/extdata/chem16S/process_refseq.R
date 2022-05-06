# chem16S/process_refseq.R

# Functions to create intermediate data files from RefSeq amino acid compositions

# Calculate amino acid composition of taxonomic groups at genus and higher ranks --> taxon_AA.csv
# mkAA()       
# Calculate chemical metrics (nH2O, ZC) for each RefSeq group --> taxon_metrics.csv
# mkmetrics()         

# Calculate average amino acid composition of genus and higher ranks in RefSeq 20200911
mkAA <- function(ranks = c("genus", "family", "order", "class", "phylum", "superkingdom")) {

  # Read RefSeq amino acid compositions and taxon names
  refseq <- read.csv(system.file("extdata/refseq/protein_refseq.csv.xz", package = "JMDplots"), as.is = TRUE)
  taxa <- read.csv(system.file("extdata/refseq/taxid_names.csv.xz", package = "JMDplots"), as.is = TRUE)
  # Make sure the data tables have consistent taxids
  stopifnot(all(refseq$organism == taxa$taxid))
  # Keep taxids classified at species level 20220104
  isspecies <- !is.na(taxa$species)
  refseq <- refseq[isspecies, ]
  taxa <- taxa[isspecies, ]
  # Exclude Bacteria and Archaea species with less than 500 sequences
  ivirus <- taxa$superkingdom == "Viruses"
  ivirus[is.na(ivirus)] <- FALSE
  ilow <- refseq$chains < 500 & !ivirus
  refseq <- refseq[!ilow, ]
  taxa <- taxa[!ilow, ]
  # Divide by number of reference sequences in each species to get mean AA composition per protein 20210605
  refseq[, 5:25] <- refseq[, 5:25] / refseq$chains

  # Make a list to hold the output
  out <- vector("list", length(ranks))
  names(out) <- ranks

  # Loop over ranks
  for(rank in ranks) {

    # Find the column corresponding to this rank
    icol <- grep(rank, colnames(taxa))

    # Get names of all taxa at this rank
    names <- na.omit(unique(taxa[, icol]))
    print(paste(rank, length(names)))

    # Create blank amino acid data frame
    # NOTE: column names are chosen to be compatible with 'protein' data in CHNOSZ
    AAtmp <- structure(list(
      protein = NA, organism = NA, ref = NA, abbrv = NA, chains = NA,
      Ala = NA, Cys = NA, Asp = NA, Glu = NA, Phe = NA,
      Gly = NA, His = NA, Ile = NA, Lys = NA, Leu = NA,
      Met = NA, Asn = NA, Pro = NA, Gln = NA, Arg = NA,
      Ser = NA, Thr = NA, Val = NA, Trp = NA, Tyr = NA), row.names = 1L, class = "data.frame")
    AA <- AAtmp[rep(1, length(names)), ]
    AA$protein <- rank
    AA$organism <- names
    rownames(AA) <- 1:length(names)

    # Loop over names
    for(i in 1:length(names)) {
      # Find all RefSeq species that have this taxon name
      taxon <- names[i]
      istax <- taxa[, icol] == taxon
      istax[is.na(istax)] <- FALSE
      if(any(istax)) {
        # Sum the number of species ("chains" column) and amino acid composition for this taxon
        sumAA <- colSums(refseq[istax, 5:25])
        # Divide by the number of species to get the mean amino acid composition for this taxon 20220107
        meanAA <- sumAA / sumAA[1]
        AA[i, 5:25] <- meanAA
        # Put the number of species into the "ref" column
        AA$ref[i] <- sum(istax)
        # Put the parent taxon into the "abbrv" column
        parent <- "Root"
        if(icol < ncol(taxa)) {
          # All taxa with the same name and rank should have the same parent, but if not, list them all
          parent <- unique(taxa[istax, icol+1])
          if(length(parent) > 1) parent <- paste(parent, collapse = ";")
        }
        AA$abbrv[i] <- parent
      }
    }

    if(rank == "genus") {
      # Add individual taxids that are used for RDP-NCBI mappings 20200922
      #addspecies <- refseq$ref %in% c("Candidatus Marinimicrobia bacterium", "Luteitalea pratensis")
      # Can't use Candidatus Marinimicrobia bacterium because it only has 4 RefSeq sequences 20220126
      addspecies <- refseq$ref %in% c("Luteitalea pratensis")
      adds <- refseq[addspecies, ]
      adds$organism <- adds$ref
      adds$ref <- 1
      adds$protein <- "species"
      AA <- rbind(adds, AA)
    }

    # Round the values
    AA[, 5:25] <- round(AA[, 5:25], 2)
    out[[rank]] <- AA

  }

  # Combine the data frames for all ranks
  out <- do.call(rbind, out)
  # Replace NA parent with ""
  out$abbrv[is.na(out$abbrv)] <- ""
  write.csv(out, "taxon_AA.csv", row.names = FALSE, quote = FALSE)
}

# Compute chemical metrics for each RefSeq group 20200927
mkmetrics <- function() {
  # Read amino acid compositions of all groups
  AA <- read.csv("taxon_AA.csv", as.is = TRUE)
  # Build output data frame; rename columns for better readability
  out <- data.frame(rank = AA$protein, group = AA$organism, ntaxa = AA$ref, parent = AA$abbrv, nH2O = NA, ZC = NA, nC = NA)
  # Calculate metrics
  out$nH2O <- round(canprot::H2OAA(AA), 6)
  out$ZC <- round(canprot::ZCAA(AA), 6)
  out$nC <- round(CAA(AA), 6)
  write.csv(out, "taxon_metrics.csv", row.names = FALSE, quote = FALSE)
}

# Function used in mkmetrics() to calculate number of carbon atoms in amino acid compositions 20200927
CAA <- function(AAcomp) {
  # the number of carbons of the amino acids
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  iAA <- match(colnames(AAcomp)[isAA], names(nC_AA))
  # calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA]) * nC_AA[iAA])
  # calculate the total nC, then the per-residue nC
  nCtot <- rowSums(multC)
  nCtot / rowSums(AAcomp[, isAA])
}

