# chem16S/RefSeq/taxon_AA.R

# Functions to generate amino acid compositions and chemical metrics of
# higher-level taxa from RefSeq species-level reference proteomes

# Notes moved from README.md on 20230708
#* Only taxids classified at the species level are used, and archaeal and bacterial species with less than 500 reference protein sequences are excluded;
#* For each species-level taxid, the total amino acid composition is converted to per-protein mean amino acid composition (this is done so that species with different proteome sizes contribute equally to the reference proteomes of higher-level taxa);
#* For each genus, the mean amino acid compositions of all species-level taxids in that genus are summed and divided by the number of taxids to get the amino acid composition of the reference proteome;
#* Analogously, the mean amino acid compositions of all species-level taxids in each family, order, class, and phylum are used to get the reference proteomes for taxa at those levels.

# Calculate amino acid composition of taxonomic groups at genus and higher ranks --> taxon_AA.csv
# taxon_AA()       

# Calculate average amino acid composition of genus and higher ranks in RefSeq 20200911
taxon_AA <- function(ranks = c("genus", "family", "order", "class", "phylum", "superkingdom")) {

  # Read RefSeq amino acid compositions and taxon names
  refseq <- read.csv(system.file("RefDB/RefSeq_206/genome_AA.csv.xz", package = "JMDplots"), as.is = TRUE)
  taxa <- read.csv(system.file("RefDB/RefSeq_206/taxonomy.csv.xz", package = "chem16S"), as.is = TRUE)
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
