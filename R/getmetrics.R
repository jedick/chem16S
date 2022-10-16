# chem16S/getmetrics.R
# Get chemical metrics from RDP classifications and RefSeq reference proteomes 20200927
# Moved to chem16S 20220505

getmetrics <- function(RDP = NULL, map = NULL, taxon_AA = NULL, groups = NULL) {

  # Exclude NA mappings
  RDP <- RDP[!is.na(map), ]
  map <- na.omit(map)
  if(length(map) == 0) stop("no available mappings RefSeq taxa!")

  # Get amino acid compositions of taxa compiled from RefSeq sequences
  # (no longer using precompiled metrics in taxon_metrics.csv 20220108)
  AAfile <- system.file("extdata/RefSeq/taxon_AA.csv", package = "chem16S")
  if(is.null(taxon_AA)) taxon_AA <- read.csv(AAfile, as.is = TRUE)
  # Keep only those taxa used in the mapping
  taxon_AA <- taxon_AA[map, ]
  equalrank <- RDP$rank == taxon_AA$protein
  # Don't test particular RDP-NCBI mappings that cross ranks
  iclassCyano <- RDP$rank == "class" & RDP$name == "Cyanobacteria"
  igenusSparto <- RDP$rank == "genus" & RDP$name == "Spartobacteria_genera_incertae_sedis"
  iclassActino <- RDP$rank == "class" & RDP$name == "Actinobacteria"
  igenusVerruco <- RDP$rank == "genus" & RDP$name == "Subdivision3_genera_incertae_sedis"
  igenusGpI <- RDP$rank == "genus" & RDP$name == "GpI"
  igenusGpIIa <- RDP$rank == "genus" & RDP$name == "GpIIa"
  equalrank <- equalrank[!(iclassCyano | igenusSparto | iclassActino | igenusVerruco | igenusGpI | igenusGpIIa)]
  stopifnot(all(equalrank))

  # Get classification matrix (rows = taxa, columns = samples)
  RDPmat <- as.matrix(RDP[, -(1:4), drop = FALSE])
  # Get amino acid matrix (rows = taxa, columns = amino acid)
  AAmat <- as.matrix(taxon_AA[, 6:25, drop = FALSE])

  if(is.null(groups)) {
    # Calculate amino acid composition for each sample (weighted by taxon abundances)
    AAcomp <- t(RDPmat) %*% AAmat
    # Calculate ZC and nH2O from amino acid composition
    ZC <- ZCAA(AAcomp)
    nH2O <- H2OAA(AAcomp)
    # Create output data frame
    out <- data.frame(Run = colnames(RDPmat), nH2O = nH2O, ZC = ZC)
  } else {
    # Split data into sample groups and calculate metrics for each group 20210607
    nH2O <- ZC <- numeric()
    for(i in 1:length(groups)) {
      # Use rowSums to combine all samples in each group into one meta-sample
      thisRDP <- rowSums(RDPmat[, groups[[i]], drop = FALSE])
      AAcomp <- t(thisRDP) %*% AAmat
      ZC <- c(ZC, ZCAA(AAcomp))
      nH2O <- c(nH2O, H2OAA(AAcomp))
    }
    sample <- names(groups)
    if(is.null(sample)) sample <- 1:length(groups)
    out <- data.frame(Run = rep(NA, length(groups)), sample = sample, nH2O = nH2O, ZC = ZC)
  }

  rownames(out) <- 1:nrow(out)
  out

}

