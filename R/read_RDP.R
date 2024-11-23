# chem16S/read_RDP.R
# Calculate chemical metrics based on 16S data and RefSeq proteins 20200902
# Revised to include "unclassified" groups in RDP (i.e. classified above genera) 20200911
# Moved to JMDplots 20210416
# Moved to chem16S 20220505

# Read and filter RDP results for all samples in a study 20200912
read_RDP <- function(file, lineage = NULL, mincount = 200, lowest.level = NULL, drop.groups = FALSE, quiet = FALSE) {

  # Read file
  dat <- read.table(file, sep = "\t", header = TRUE, check.names = FALSE)

  # Basename of file to use in messages 20220505
  basefile <- basename(file)
  basetxt <- paste0(" [", basefile, "]")

  # Columns with classification counts
  icol <- 5:ncol(dat)

  # Rows with "unclassified" assignments
  iunclass <- grepl("unclassified_", dat$name)
  if(any(iunclass)) {
    # Find the ranks for the "unclassified" sequences
    # (which means classified at higher rank than genus)
    unclassname <- gsub("unclassified_", "", dat$name[iunclass])
    # Split the lineage text and find the length of each one
    slineage <- strsplit(dat$lineage[iunclass], ";")
    sllength <- vapply(slineage, length, 1)
    # The rank is a certain number of positions before the end
    irank <- sllength - 2
    unclassrank <- mapply("[", slineage, irank)
    dat$rank[iunclass] <- unclassrank
    dat$name[iunclass] <- unclassname
  }
  if(dat$lineage[1] == "null") {
    # This is for the actual output from the RDP Classifier
    # We want to keep genus-level classifications; also:
    # Remove the counts of classifications at higher ranks (which are also counted at lower ranks)
    # *except* for the ones labeled "unclassified" (which are not counted at lower ranks)
    igenus <- dat$rank == "genus"
    out <- dat[igenus | iunclass, ]
    isRDP <- TRUE
  } else {
    # Otherwise, the taxon counts are assembled from SI or OTU tables; use all the counts
    out <- dat
    isRDP <- FALSE
  }
  # Get the total counts
  totalcounts <- colSums(out[, icol, drop = FALSE])
  if(isRDP) {
    # Get the "rootrank" counts
    rootcounts <- dat[1, icol]
    # Stop if total counts doesn't equal "rootrank" counts
    stopifnot(all(abs(totalcounts - rootcounts) < 0.1))
  }

  # Keep specified lineage 20200924
  if(!is.null(lineage)) {
    precount <- sum(out[, icol])
    irow <- grepl(lineage, out$lineage)
    if(!any(irow)) stop(paste("nothing available for lineage =", lineage))
    out <- out[irow, , drop = FALSE]
    if(!quiet) {
      postcount <- sum(out[, icol])
      lpercent <- formatC(postcount / precount * 100, 1, format = "f")
      print(paste0("read_RDP", basetxt, ": keeping ", lineage, " lineage (", lpercent, "%)"))
    }
    # Recalculate total counts 20211008
    totalcounts <- colSums(out[, icol, drop = FALSE])
  }

  # Keep the rows with any counts > 0
  groupcounts <- rowSums(out[, icol, drop = FALSE])
  out <- out[groupcounts > 0, , drop = FALSE]

  if(!quiet) {
    # Get the number of counts classified at genus level
    igenus <- out$rank == "genus"
    genuscounts <- colSums(out[igenus, icol, drop = FALSE], na.rm = TRUE)
    # Print percentage of assignments at genus level
    genuspercent <- round(100 * sum(genuscounts) / sum(totalcounts))
    print(paste0("read_RDP", basetxt, ": ", genuspercent, "% of classifications at genus level"))
  }

  # Drop groups above phylum or not archaea or bacteria
  if(drop.groups) {
    RDPgroups <- paste(out$rank, out$name, sep = "_")
    RDPgroups[grepl("Eukaryota", out$lineage)] <- "Eukaryota"
    rmgroups <- c(
      # Remove classifications at root and domain level (Bacteria and Archaea) 20200922
      "rootrank_Root", "domain_Bacteria", "domain_Archaea",
      # Chloroplast 20200922
      "class_Chloroplast", "family_Chloroplast", "genus_Chlorophyta", "genus_Bacillariophyta",
      # Eukaryota 20211012
      "Eukaryota"
    )
    isrm <- RDPgroups %in% rmgroups
    if(any(isrm)) {
      irm <- which(isrm)
      rmpercent <- round(rowSums(out[irm, icol, drop = FALSE]) / sum(totalcounts) * 100, 1)
      for(i in seq_along(irm)) {
        # Only print message if removed group is >= 0.1% 20200927
        if(!quiet & rmpercent[i] >= 0.1) print(paste0("read_RDP", basetxt, ": removing ", RDPgroups[irm[i]], " (", rmpercent[i], "%)"))
      }
      out <- out[!isrm, , drop = FALSE] 
    }
    # Recalculate total counts
    totalcounts <- colSums(out[, icol, drop = FALSE])
  }

  # Discard samples with < mincount total counts 20201001
  ismall <- totalcounts < mincount
  if(any(ismall)) {
    if(!quiet) print(paste0("read_RDP", basetxt, ": discarding ", sum(ismall), " samples with < ", mincount, " total counts"))
    out <- out[, c(TRUE, TRUE, TRUE, TRUE, !ismall)]
    if(ncol(out) == 4) stop("No samples are left! (try lowering 'mincount')")
    icol <- 5:ncol(out)
  }

  if(!quiet) {
    # Recalculate total counts
    totalcounts <- colSums(out[, icol, drop = FALSE])
    # Report the median number of counts 20200917
    # Change this to range 20200924
    if(length(totalcounts) == 0) {
      print(paste0("read_RDP", basetxt, ": no samples contain at least ", mincount, " counts"))
    } else {
      print(paste0("read_RDP", basetxt, ": range of counts is ", paste(round(range(totalcounts)), collapse = " to "), ""))
    }
  }

#  # Adjust counts for 16S rRNA gene copy number 20200927
#  if(cn) {
#    # Data from rdp_classifier_2.13/src/data/classifier/16srrna/bergeyTrainingTree.xml (20200720)
#    bergey <- read.csv(system.file("extdata/RDP/bergeyTrainingTree.csv", package = "chem16S"))
#    # Paste together rank and name
#    RDPgroups <- paste(out$rank, out$name, sep = "_")
#    bergeygroups <- paste(bergey$rank, bergey$name, sep = "_")
#    # Get the copy number for these groups
#    ibergey <- match(RDPgroups, bergeygroups)
#    cpNumber <- bergey$cpNumber[ibergey]
#    # Divide RDP Classifier counts by copy number
#    # Round to 3 decimal places following output of RDP Classifier with copy-number adjustment
#    out[, icol] <- round(out[, icol] / cpNumber, 3)
#  }

  RDP <- out

  # Trim classifications to specified lowest level 20220112
  if(!is.null(lowest.level)) {
    lineage <- RDP$lineage
    # Split the lineage at the specified level
    trimmed.lineage <- sapply(strsplit(RDP$lineage, lowest.level), "[", 1)
    # Get the name of the taxon at this level
    new.name <- sapply(strsplit(trimmed.lineage, ";"), tail, 1)
    # To make the new lineage, add the name of the rank at the end
    new.lineage <- paste0(trimmed.lineage, lowest.level, ";")
    # Only update lineages that have the lowest-level classification
    has.ll <- grepl(lowest.level, lineage)
    RDP$lineage[has.ll] <- new.lineage[has.ll]
    RDP$name[has.ll] <- new.name[has.ll]
    RDP$rank[has.ll] <- lowest.level
    # Check that each unique lineage string corresponds to a unique rank-name combination
    stopifnot(identical(duplicated(RDP$lineage), duplicated(paste(RDP$rank, RDP$name))))
    # Combine duplicated lineages
    newRDP <- aggregate(RDP[, icol], list(lineage = RDP$lineage), sum)
    # Prepend taxid, lineage, name, rank columns
    iRDP <- match(newRDP$lineage, RDP$lineage)
    newRDP <- cbind(RDP[iRDP, 1:4], newRDP[, -1, drop = FALSE])
    RDP <- newRDP
  }

  # Set attributes 20220505
  # Original attributes (names, row.names, class)
  attr.orig <- attributes(RDP)
  # New attributes
  # (NOTE: NULL values are not stored)
  attr.new <- list(file = file, lineage = lineage, mincount = mincount, lowest.level = lowest.level)
  attributes(RDP) <- c(attr.orig, attr.new)

  RDP

}
