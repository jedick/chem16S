# chem16S/get_metrics.R
# Get chemical metrics from taxonomic classifications and RefSeq reference proteomes 20200927
# Moved to chem16S 20220505
# Add refdb argument 20221016

get_metrics <- function(RDP = NULL, map = NULL, refdb = "GTDB_220", taxon_AA = NULL,
  groups = NULL, return_AA = FALSE, zero_AA = NULL, metrics = c("Zc", "nO2", "nH2O")) {

  # Exclude NA mappings
  RDP <- RDP[!is.na(map), ]
  map <- na.omit(map)
  if(length(map) == 0) stop(paste("no available mappings to taxa in", refdb, "reference database"))

  # Get amino acid compositions of taxa compiled from RefSeq or GTDB
  AApath <- file.path("RefDB", refdb, "taxon_AA.csv.xz")
  AAfile <- system.file(AApath, package = "chem16S")
  if(is.null(taxon_AA)) taxon_AA <- read.csv(AAfile, as.is = TRUE)
  
  # Apply the mapping
  taxon_AA <- taxon_AA[map, ]

  # Get classification matrix (rows = taxa, columns = samples)
  RDPmat <- as.matrix(RDP[, -(1:4), drop = FALSE])
  # Get amino acid matrix (rows = taxa, columns = 'chains' (number of proteins) and one column for each amino acid)
  AAmat <- as.matrix(taxon_AA[, 5:25, drop = FALSE])

  if(is.null(groups)) {
    # Calculate amino acid composition for each sample (weighted by taxon abundances)
    AAcomp <- as.data.frame(t(RDPmat) %*% AAmat)
    # Check that the number of protein sequences is equal to the number of classifications
    stopifnot(all(AAcomp[, "chains"] == colSums(RDPmat)))
    # Names of samples are sequencing run names
    samplecols <- data.frame(Run = colnames(RDPmat))
  } else {
    # Aggregate samples to calculate metrics for each group 20210607
    AAcomp <- lapply(groups, function(group) {
      # Use rowSums to combine all samples in each group into one meta-sample
      thisRDP <- rowSums(RDPmat[, group, drop = FALSE])
      t(thisRDP) %*% AAmat
    })
    AAcomp <- as.data.frame(do.call(rbind, AAcomp))
    # Check that the total number of protein sequences is equal to the total number of classifications
    stopifnot(sum(AAcomp[, "chains"]) == sum(RDPmat))
    # Names of samples are group names
    group <- names(groups)
    if(is.null(group)) group <- 1:length(groups)
    samplecols <- data.frame(group = group)
  }

  # Set counts of specified amino acids to zero 20221018
  if(!is.null(zero_AA)) AAcomp[, zero_AA] <- 0
  if(return_AA) {
    # Just return the amino acid composition
    out <- cbind(samplecols, AAcomp)
  } else {
    # Divide amino acid composition by number of proteins in order to calculate
    # length and other protein-level metrics (pMW, etc.) correctly 20240329
    AAcomp <- AAcomp / AAcomp$chains
    # Calculate chemical metrics from amino acid composition
    metrics_values <- calc_metrics(AAcomp, metrics)
    colnames(metrics_values) <- metrics
    # Create output data frame
    out <- cbind(samplecols, metrics_values)
  }

  rownames(out) <- 1:nrow(out)
  out

}
