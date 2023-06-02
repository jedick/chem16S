# chem16S/getmetrics.R
# Get chemical metrics from RDP classifications and RefSeq reference proteomes 20200927
# Moved to chem16S 20220505
# Add refdb argument 20221016

getmetrics <- function(RDP = NULL, map = NULL, refdb = "RefSeq", taxon_AA = NULL, groups = NULL, return_AA = FALSE, zero_AA = NULL) {

  # Exclude NA mappings
  RDP <- RDP[!is.na(map), ]
  map <- na.omit(map)
  if(length(map) == 0) stop(paste("no available mappings to taxa in", refdb, "reference database"))

  # Get amino acid compositions of taxa compiled from:
  #   - RefSeq (no longer using precompiled metrics in taxon_metrics.csv 20220108) or
  #   - GTDB 20221016
  AApath <- file.path("extdata", refdb, "taxon_AA.csv.xz")
  AAfile <- system.file(AApath, package = "chem16S")
  if(is.null(taxon_AA)) taxon_AA <- read.csv(AAfile, as.is = TRUE)
  
  # Apply the mapping
  taxon_AA <- taxon_AA[map, ]
  # Check that the before and after ranks are equal
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
  AAmat <- as.matrix(taxon_AA[, 5:25, drop = FALSE])

  if(is.null(groups)) {
    # Calculate amino acid composition for each sample (weighted by taxon abundances)
    AAcomp <- t(RDPmat) %*% AAmat
    # Check that the number of protein sequences is equal to the number of classifications
    stopifnot(all(AAcomp[, "chains"] == colSums(RDPmat)))
    # Set counts of specified amino acids to zero 20221018
    if(!is.null(zero_AA)) AAcomp[, zero_AA] <- 0
    # Calculate ZC, nO2, and nH2O from amino acid composition
    ZC <- ZCAA(AAcomp)
    nO2 <- O2AA(AAcomp)
    nH2O <- H2OAA(AAcomp)
    # Create output data frame
    out <- data.frame(Run = colnames(RDPmat), ZC = ZC, nO2 = nO2, nH2O = nH2O)
    if(return_AA) out <- cbind(data.frame(Run = colnames(RDPmat)), AAcomp)
  } else {
    # Split data into sample groups and calculate metrics for each group 20210607
    nH2O <- nO2 <- ZC <- numeric()
    for(i in 1:length(groups)) {
      # Use rowSums to combine all samples in each group into one meta-sample
      thisRDP <- rowSums(RDPmat[, groups[[i]], drop = FALSE])
      AAcomp <- t(thisRDP) %*% AAmat
      # Set counts of specified amino acids to zero 20221018
      if(!is.null(zero_AA)) AAcomp[, zero_AA] <- 0
      ZC <- c(ZC, ZCAA(AAcomp))
      nO2 <- c(nO2, O2AA(AAcomp))
      nH2O <- c(nH2O, H2OAA(AAcomp))
    }
    group <- names(groups)
    if(is.null(group)) group <- 1:length(groups)
    out <- data.frame(group = group, ZC = ZC, nO2 = nO2, nH2O = nH2O)
    if(return_AA) out <- cbind(data.frame(group = group), AAcomp)
  }

  rownames(out) <- 1:nrow(out)
  out

}
