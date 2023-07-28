# chem16S/get_metadata.R
# Get metadata for a study, appending columns for pch and col 20200914
# Moved to chem16S 20220505

get_metadata <- function(file, metrics = NULL) {
  # Read metadata file
  metadata <- read.csv(file, as.is = TRUE, check.names = FALSE)

  # Use NULL pch as flag for unavailable dataset 20210820
  pch <- NULL
  # Basename of file to use for identifying study 20220505
  basefile <- basename(file)

  # Assign pch and col based on specific metadata for each study

  if(basefile == "BGPF13.csv") {
    # Bowen De León et al., (2013)
    # Heart Lake Geyser Basin, Yellowstone
    pch <- sapply(metadata$domain, switch, Bacteria = 22, Archaea = 23)
    col <- sapply(metadata$domain, switch, Bacteria = 5, Archaea = 6)
  }
  if(basefile == "SMS+12.csv") {
    # Swingley et al. (2012)
    # Bison Pool, Yellowstone
    pch <- ifelse(metadata$O2 > 0.5, 24, 25)
    col <- ifelse(metadata$O2 > 0.5, 4, 2)
  }
  if(basefile == "HLA+16.csv") {
    # Herlemann et al. (2016)
    # Baltic Sea
    type <- rep("moderate", nrow(metadata))
    type[metadata$salinity < 6] <- "low"
    type[metadata$salinity > 20] <- "high"
    pch <- sapply(type, switch, low = 24, moderate = 20, high = 21)
    col <- sapply(type, switch, low = 3, moderate = 1, high = 4)
  }

  if(is.null(pch)) stop(paste(basefile, "exists, but is not set up for processing"))

  metadata <- cbind(metadata, pch, col)

  # Return both metadata and metrics, if provided 20220506
  if(is.null(metrics)) metadata else {
    # Keep metadata only for samples with metrics 20201006
    metadata <- metadata[metadata$Run %in% metrics$Run, ]
    # Put metrics in same order as metadata 20220505
    imet <- match(metadata$Run, metrics$Run)
    metrics <- metrics[imet, ]
    # Insert sample column in metrics
    # Use first column name starting with "sample" or "Sample" 20210818
    sampcol <- grep("^sample", colnames(metadata), ignore.case = TRUE)[1]
    metrics <- data.frame(Run = metrics[, 1], sample = metadata[, sampcol], metrics[, -1])
    list(metadata = metadata, metrics = metrics)
  }

}
