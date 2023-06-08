# chem16S/physeq.R
# Calculating chemical metrics from phyloseq objects

# Returns data frame with lowest-level classifications for each OTU 20230607
# (genus to domain level - use column names similar to output from RDP Classifier)
physeq_taxacounts <- function(physeq, split = TRUE) {

  # Get taxonomy table
  tax <- tax_table(physeq)

  # Get OTU table
  # NOTE: We want to keep taxa as rows (the opposite of estimate_richness())
  if(!split) {
    # Sum the taxonomic abundances
    OTU <- as.data.frame(taxa_sums(physeq))
  } else {
    OTU <- otu_table(physeq)
    if(!taxa_are_rows(physeq)) OTU <- t(OTU)
  }

  # Convert OTU table to data frame
  taxacounts <- as.data.frame(OTU)
  # Initialize taxid (OTU name), lineage, name, and rank columns
  taxid <- rownames(taxacounts)
  taxacounts <- cbind(taxid, lineage = NA, name = NA, rank = NA, taxacounts)
  # Loop over taxonomic ranks
  for(rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
    # Get classifications at this rank
    names <- as.vector(tax[, rank])
    # Insert non-NA classifications into data frame
    is.classified <- !is.na(names)
    taxacounts$name[is.classified] <- names[is.classified]
    taxacounts$rank[is.classified] <- tolower(rank)
    # Add name to lineage
    has.lineage <- !is.na(taxacounts$lineage)
    # If lineage is NA, just use the first available classification (e.g. "Archaea")
    taxacounts$lineage[is.classified & !has.lineage] <- names[is.classified]
    # If lineage is not NA, append current classification to existing lineage
    ichl <- is.classified & has.lineage
    taxacounts$lineage[ichl] <- paste(taxacounts$lineage[ichl], names[ichl], sep = ";")
  }
  taxacounts

}

# Calculate chemical metrics from phyloseq object
physeq_metrics <- function(physeq, split = TRUE, metrics = c("ZC", "nO2", "nH2O"), quiet = TRUE) {
  # Obtain data frame with lowest-level (to genus) classifications for each OTU
  taxacounts <- physeq_taxacounts(physeq, split = split)
  # Map names to NCBI taxonomy
  map <- mapRDP(taxacounts, quiet = quiet)
  # Calculate chemical metrics for community reference proteomes
  met <- getmetrics(taxacounts, map)
  # Put sample names in rownames (analogous to phyloseq::estimate_richness)
  rownames(met) <- met$Run
  met <- met[, -1]
  # Keep only selected metrics
  met <- met[, metrics, drop = FALSE]
  met
}

plot_metrics <- function(physeq, x = "samples", color = NULL, shape = NULL, title = NULL,
  scales = "free_y", nrow = 1, metrics = c("ZC", "nO2", "nH2O"), sortby = NULL) { 

  # Calculate the chemical metrics
  pmDF <- physeq_metrics(physeq, split = TRUE, metrics = metrics)

  # Make the plotting data.frame.
  # This coerces to data.frame, required for reliable output from reshape2::melt()
  if( !is.null(sample_data(physeq, errorIfNULL = FALSE)) ){
    # Include the sample data, if it is there.
    DF <- data.frame(pmDF, sample_data(physeq))
  } else {
    # If no sample data, leave it out.
    DF <- data.frame(pmDF)
  }

  if( !"samples" %in% colnames(DF) ){
    # If there is no "samples" variable in DF, add it
    DF$samples <- sample_names(physeq)
  }
  # sample_names used to be default, and should also work.
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
  mdf = melt(DF, metric.vars = metrics)

#  # Initialize the se column. Helpful even if not used.
#  mdf$se <- NA_integer_

  ## Interpret metrics
  # If not provided (default), keep all 
  if( !is.null(metrics) ){
    if( any(metrics %in% as.character(mdf$variable)) ){
      # If any metrics were in mdf, then subset to just those.
      mdf <- mdf[as.character(mdf$variable) %in% metrics, ]
    } else {
      # Else, print warning about bad option choice for metrics, keeping all.
      warning("Argument to `metrics` not supported. All chemical metrics (should be) included in plot.")
    }
  }

  # Address `sortby` argument
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
  #metrics_map <- aes_string(x = x, y = "value", colour = color, shape = shape)

  # aes_string is deprecated - use tidy evaluation instead 20230608 jmd
  # TODO: Find elegant way to handle NULL values ...
  #       See https://stackoverflow.com/questions/57350795/r-rlang-handle-null-arguments
  if(!is.null(color) & !is.null(shape)) {
    metrics_map <- aes(.data[[x]], .data[["value"]], colour = .data[[color]], shape = .data[[shape]])
  } else if(!is.null(color)) {
    metrics_map <- aes(.data[[x]], .data[["value"]], colour = .data[[color]])
  } else if(!is.null(shape)) {
    metrics_map <- aes(.data[[x]], .data[["value"]], shape = .data[[shape]])
  } else {
    metrics_map <- aes(.data[[x]], .data[["value"]])
  }

  # Make the ggplot.
  p <- ggplot(mdf, metrics_map) + geom_point(na.rm = TRUE)

#  # Add error bars if mdf$se is not all NA
#  if( any(!is.na(mdf[, "se"])) ){
#    p <- p + geom_errorbar(aes(ymax = value + se, ymin = value - se), width = 0.1) 
#  }

  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))
  # Add y-label 
  p <- p + ylab("Chemical Metric") 
  # Facet wrap using user-options
  p <- p + facet_wrap(~variable, nrow = nrow, scales = scales)
  # Optionally add a title to the plot
  if( !is.null(title) ){
    p <- p + ggtitle(title)
  }
  return(p)

}
