# chem16S/get_metric_byrank.R
# Get a single chemical metric for taxa in each sample aggregated to specified rank
# First version adapted from get_metrics() 20241119

get_metric_byrank <- function(RDP = NULL, map = NULL, refdb = "GTDB_220", taxon_AA = NULL,
  groups = NULL, zero_AA = NULL, metric = "Zc", rank = "genus") {

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

  # Get lineages of all taxa into table format 20241119
  lineage_table <- t(sapply(strsplit(RDP$lineage, ";", fixed = TRUE), "[", seq(1, 13, 2)))
  ranks <- c("rootrank", "domain", "phylum", "class", "order", "family", "genus")
  lineage_table <- data.frame(lineage_table)
  colnames(lineage_table) <- ranks
  # Fill in NA names with higher-level names
  for(i in 2:length(ranks)) {
    ina <- is.na(lineage_table[, i])
    lineage_table[ina, i] <- paste(lineage_table[ina, i-1], ranks[i], sep = "_")
  }

  if(!is.null(groups)) {
    # Aggregate samples to calculate metrics for each group 20210607
    RDPmat <- sapply(groups, function(group) rowSums(RDPmat[, group, drop = FALSE]))
    # Add group names if missing
    if(is.null(colnames(RDPmat))) colnames(RDPmat) <- 1:ncol(RDPmat)
  }

  # For each sample, weight amino acid composition by
  # taxon abundance and aggregate to specified rank 20241119
  taxa_metric_list <- lapply(1:ncol(RDPmat), function(isample) {

    # Multiply amino acid composition by taxon counts
    weighted_AA <- AAmat * RDPmat[, isample]
    # Prepend taxon names starting with lowest level (genus)
    AAcomp <- cbind(data.frame(taxon = lineage_table$genus), weighted_AA)

    # Number of target rank (1 - genus, 7 - rootrank)
    nrank <- 8 - match(rank, ranks)
    # Aggregate amino acid compositions to higher ranks
    for(irank in 1:nrank) {
      # Do nothing for genus
      if(irank == 1) next
      # Get the indices of the lower rank
      ilineage <- match(AAcomp$taxon, lineage_table[, 9 - irank])
      if(any(is.na(ilineage))) stop(paste("missing name in lineage table:", AAcomp$taxon[is.na(ilineage)]))
      # Get the names of the higher rank
      higher_name <- lineage_table[ilineage, 8 - irank]
      # Update the names in the dataframe
      AAcomp$taxon <- higher_name
      # Aggregate to this rank
      AAcomp <- aggregate(. ~ taxon, AAcomp, sum)
    }

    # Set counts of specified amino acids to zero 20221018
    if(!is.null(zero_AA)) AAcomp[, zero_AA] <- 0
    # Divide amino acid composition by number of proteins in order to calculate
    # length and other protein-level metrics (pMW, etc.) correctly 20240329
    AAcomp[, -1] <- AAcomp[, -1] / AAcomp$chains
    # Calculate chemical metric from amino acid composition
    taxa_metric <- calc_metrics(AAcomp, metric)
    rownames(taxa_metric) <- AAcomp$taxon
    taxa_metric

  })

  # Convert list to dataframe and add sample names
  taxa_metric <- do.call(cbind, taxa_metric_list)
  colnames(taxa_metric) <- colnames(RDPmat)
  # Return transpose (samples on rows)
  data.frame(t(taxa_metric))

}
