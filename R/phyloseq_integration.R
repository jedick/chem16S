# chem16S/phyloseq_integration.R
# Calculating chemical metrics from phyloseq objects

## Required packages:
#library(phyloseq)
#library(reshape2)
#library(ggplot2)
#library(plyr)

# Returns data frame with lowest-level classifications for each OTU 20230607
# (genus to domain level - use column names similar to output from RDP Classifier)
ps_taxacounts <- function(physeq, split = TRUE) {

  # Get taxonomy table
  taxtable <- phyloseq::tax_table(physeq)

  # Get OTU table
  # NOTE: Here we put taxa in rows (the opposite of phyloseq::estimate_richness())
  if(!split) {
    # Sum the taxonomic abundances
    OTU <- phyloseq::taxa_sums(physeq)
  } else {
    OTU <- phyloseq::otu_table(physeq)
    if(!phyloseq::taxa_are_rows(physeq)) OTU <- t(OTU)
  }

  # Convert OTU table to data frame
  taxacounts <- as.data.frame(OTU, check.names = FALSE)
  # Initialize taxid (OTU name), lineage, name, and rank columns
  taxid <- rownames(taxacounts)
  taxacounts <- cbind(taxid, lineage = NA, name = NA, rank = NA, taxacounts)
  # Get ranks from column names of taxonomic table 20240121
  ps_ranks <- phyloseq::rank_names(physeq)
  # Allow Domain or Kingdom, and exclude Species
  allowed_ranks <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  ranks <- intersect(ps_ranks, allowed_ranks)
  # Loop over taxonomic ranks
  for(rank in ranks) {
    # Get classifications for this rank
    names <- as.vector(taxtable[, rank])
    # Insert non-NA classifications into data frame
    is.classified <- !is.na(names)
    taxacounts$name[is.classified] <- names[is.classified]
    taxacounts$rank[is.classified] <- tolower(rank)
    # Add name to lineage
    # If lineage is NA, start with the first available classification (e.g. "Archaea")
    has.lineage <- !is.na(taxacounts$lineage)
    taxacounts$lineage[is.classified & !has.lineage] <- names[is.classified]
    # If lineage is not NA, append current classification to existing lineage
    ichl <- is.classified & has.lineage
    taxacounts$lineage[ichl] <- paste(taxacounts$lineage[ichl], names[ichl], sep = ";")
  }

  # Exclude domain/kingdom-level classifications 20230618
  taxacounts <- taxacounts[! taxacounts$rank %in% c("domain", "kingdom"), ]
  # Exclude NA classifications 20230625
  taxacounts <- taxacounts[!is.na(taxacounts$lineage), ]

  # Sum counts for each unique lineage 20230614
  idup <- duplicated(taxacounts$lineage)
  dedup <- taxacounts[!idup, ]
  for(i in 1:nrow(dedup)) {
    lineage <- dedup$lineage[i]
    ilin <- taxacounts$lineage == lineage
    if(sum(ilin) > 1) {
      # For lineages that appear more than once, sum the counts for all occurrences
      sumcounts <- colSums(taxacounts[ilin, -(1:4)])
      dedup[i, -(1:4)] <- sumcounts
      # Paste together the taxids (sample names)
      dedup$taxid[i] <- paste(taxacounts$taxid[ilin], collapse = ";")
    }
  }

  dedup

}

# Calculate chemical metrics from phyloseq object
ps_metrics <- function(physeq, split = TRUE, refdb = "GTDB_220", quiet = FALSE, ...) {
  # Obtain data frame with lowest-level (to genus) classifications for each OTU
  taxacounts <- ps_taxacounts(physeq, split = split)
  # Map names to NCBI taxonomy
  map <- map_taxa(taxacounts, refdb = refdb, quiet = quiet)
  # Calculate chemical metrics for community reference proteomes
  met <- get_metrics(taxacounts, map, refdb = refdb, ...)
  moreargs <- list(...)
  if(isTRUE(moreargs$return_AA)) return(met)
  # Put samples in rows (analogous to phyloseq::estimate_richness())
  rownames(met) <- met$Run
  met[, -1, drop = FALSE]
}

# Plot individual chemical metrics 20230608
# Adapted by Jeffrey Dick from phyloseq::plot_richness() by Paul J. McMurdie
plot_ps_metrics <- function(physeq, metrics = c("Zc", "nO2", "nH2O"), x = "samples",
  color = NULL, shape = NULL, title = NULL,
  scales = "free_y", nrow = 1, sortby = NULL, ...) { 

  # Calculate the chemical metrics
  #pmDF <- ps_metrics(physeq, metrics = metrics, refdb = refdb, quiet = quiet)
  pmDF <- ps_metrics(physeq, metrics = metrics, ...)

  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE)) ){
    # Include the sample data, if it is there.
    DF <- data.frame(pmDF, phyloseq::sample_data(physeq), check.names = FALSE)
  } else {
    # If no sample data, leave it out
    DF <- data.frame(pmDF, check.names = FALSE)
  }

  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- phyloseq::sample_names(physeq)
  }
  # sample_names used to be default, and should also work
  # #backwardcompatibility
  if( !is.null(x) ){
    if( x %in% c("sample", "samples", "sample_names", "sample.names") ){
      x <- "samples"
    }
  } else {
    # If x was NULL for some reason, set it to "samples"
    x <- "samples"
  }
  # melt to display different chemical metrics separately
  mdf = melt(DF, measure.vars = metrics)

  # Handle `sortby` argument
  if(!is.null(sortby)){
    if(!all(sortby %in% levels(mdf$variable))){
      warning("`sortby` argument not among `metrics`. Ignored.")
    }
    if(!is.discrete(mdf[, x])){
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if(all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, x])){
      # Replace x-factor with same factor that has levels re-ordered according to `sortby`
      wh.sortby <- which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x],
                         levels = names(sort(tapply(X = mdf[wh.sortby, "value"],
                                                    INDEX = mdf[wh.sortby, x],
                                                    mean,
                                                    na.rm = TRUE, simplify = TRUE))))
    }
  }

  # Define variable mapping
  #metrics_map <- aes_string(x = x, y = "value", color = color, shape = shape)
  # aes_string is deprecated - use tidy evaluation instead 20230608 jmd
  # Start with just x and y in case color and/or shape are NULL 20230709
  metrics_map <- aes(.data[[x]], .data[["value"]])
  # https://stackoverflow.com/questions/20084104/combine-merge-two-ggplot-aesthetics
  if(!is.null(color)) metrics_map <- modifyList(metrics_map, aes(color = .data[[color]]))
  if(!is.null(shape)) metrics_map <- modifyList(metrics_map, aes(shape = .data[[shape]]))

  # Make the ggplot
  p <- ggplot(mdf, metrics_map) + geom_point(na.rm = TRUE)

  if(length(metrics) == 1) {
    # Don't use facets for a single metric 20230707
    # Label the y-axis with the formatted label
    p <- p + ylab(chemlab(metrics))
  } else {
    ## Rotate horizontal axis labels, and adjust
    #p <- p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))
    # Add y-label 
    p <- p + ylab("Chemical Metric") 
    # Define function for facet labels 20230617
    # https://stackoverflow.com/questions/3472980/how-to-change-facet-labels [outdated]
    metrics_labels <- function(variable) {
      lapply(variable, chemlab)
    }
    # Facet wrap using user-options
    # Plot expressions with label_parsed 20230617 jmd
    # https://stackoverflow.com/questions/37089052/r-ggplot2-facet-grid-how-include-math-expressions-in-few-not-all-labels
    p <- p + facet_wrap(. ~ variable, nrow = nrow, scales = scales, labeller = labeller(.default = label_parsed, variable = metrics_labels))
  }

  # Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }

  return(p)

}

# Plot two chemical metrics against each other 20230617
# Parts of this function were adapted by Jeffrey Dick from phyloseq::plot_richness() by Paul J. McMurdie
plot_ps_metrics2 <- function(physeq, metrics = c("Zc", "nH2O"), color = NULL, shape = NULL,
  title = NULL, refdb = "GTDB_220", quiet = FALSE) { 

  if(length(metrics) < 2) stop("please supply the names of two metrics in 'metrics'")
  if(length(metrics) > 2) warning("plot_ps_metrics2: 'metrics' has length > 2; using the first two")

  # Calculate the chemical metrics
  pmDF <- ps_metrics(physeq, metrics = metrics, refdb = refdb, quiet = quiet)

  # Make the plotting data.frame
  if( !is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE)) ){
    # Include the sample data, if it is there
    DF <- data.frame(pmDF, phyloseq::sample_data(physeq), check.names = FALSE)
  } else {
    # If no sample data, leave it out
    DF <- data.frame(pmDF, check.names = FALSE)
  }

  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- phyloseq::sample_names(physeq)
  }

  # Start with just x and y variables in case color and/or shape are NULL 20230709
  metrics_map <- aes(.data[[metrics[1]]], .data[[metrics[2]]])
  # https://stackoverflow.com/questions/20084104/combine-merge-two-ggplot-aesthetics
  if(!is.null(color)) metrics_map <- modifyList(metrics_map, aes(color = .data[[color]]))
  if(!is.null(shape)) metrics_map <- modifyList(metrics_map, aes(shape = .data[[shape]]))

  # Make the ggplot
  p <- ggplot(DF, metrics_map) + geom_point(na.rm = TRUE) + xlab(chemlab(metrics[1])) + ylab(chemlab(metrics[2]))

  # Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }

  return(p)

}
