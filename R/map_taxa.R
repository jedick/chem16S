# chem16S/map_taxa.R
# Map RDP to RefSeq taxonomy 20200912
# Moved to chem16S 20220505
# Add refdb argument to use RefSeq or GTDB 20221016
# TODO: warn when mapped percentage is below a certain value 20230615

map_taxa <- function(taxacounts = NULL, refdb = "GTDB_214", taxon_AA = NULL, quiet = FALSE) {

  # Make group names by combining rank and name
  INPUTgroups <- paste(taxacounts$rank, taxacounts$name, sep = "_")
  # Calculate group abundances for displaying messages
  groupcounts <- as.numeric(rowSums(taxacounts[, -(1:4), drop = FALSE]))
  # Basename of file to use in messages 20220505
  basetxt <- ""
  fileattr <- attr(taxacounts, "file")
  if(!is.null(fileattr)) basetxt <- paste0(" [", basename(fileattr), "]")

  # Manual mappings for RefSeq (NCBI) taxonomy but not for GTDB 20221016
  if(grepl("RefSeq", refdb)) {

    mapping <- getOption("manual_mappings")
    INPUT <- paste(mapping$RDP.rank, mapping$RDP.name, sep = "_")
    NCBI <- paste(mapping$NCBI.rank, mapping$NCBI.name, sep = "_")
    iswitch <- INPUTgroups %in% INPUT

    if(any(iswitch)) {
      # Make the switch
      imap <- match(INPUTgroups[iswitch], INPUT)
      INPUTgroups[iswitch] <- NCBI[imap]
      if(!quiet) {
        # Print message(s) about switched names and abundance
        from <- INPUT[imap]
        to <- NCBI[imap]
        switchcounts <- groupcounts[iswitch]
        switchpercent <- round(switchcounts / sum(groupcounts) * 100, 1)
        # Only print message for mappings of groups at least 0.1% abundant 20200927
        if(any(switchpercent >= 0.1)) {
          print(paste0("map_taxa", basetxt, ": using these manual mapping(s) to NCBI RefSeq:"))
          for(i in seq_along(from)) {
            if(switchpercent[i] >= 0.1) message(paste0(from[i], " --> ", to[i], " (", switchpercent[i], "%)"))
          }
        }
      }
    }

  }

  # Read amino acid compositions of reference proteomes for genus and higher taxonomic groups:
  #   - from RefSeq (NCBI) or
  #   - from GTDB 20221016
  AApath <- file.path("RefDB", refdb, "taxon_AA.csv.xz")
  AAfile <- system.file(AApath, package = "chem16S")
  if(is.null(taxon_AA)) {
    if(!file.exists(AAfile)) {
      available_RefDB <- dir(system.file("RefDB", package = "chem16S"))
      stop(paste0("Chosen refdb (", refdb, ") is not one of ", paste(available_RefDB, collapse = ", ")))
    }
    taxon_AA <- read.csv(AAfile, as.is = TRUE)
  }
  AAgroups <- paste(taxon_AA$protein, taxon_AA$organism, sep = "_")
  # Do the mapping!
  iAA <- match(INPUTgroups, AAgroups)
  # Get percentages of unmapped groups
  naAA <- is.na(iAA)
  nagroups <- INPUTgroups[naAA]
  nacounts <- groupcounts[naAA]
  napercent <- nacounts / sum(groupcounts) * 100
  if(!quiet) {
    # Print summary of missing groups
    naorder <- order(napercent, decreasing = TRUE)
    ordergroups <- nagroups[naorder]
    orderpercent <- round(napercent[naorder], 2)
    if(sum(naAA) > 0) namsg <- paste0("map_taxa", basetxt, ": can't map group ", ordergroups[1], " (", orderpercent[1], "%)")
    if(sum(naAA) > 1) namsg <- paste0("map_taxa", basetxt, ": can't map groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%)")
    if(sum(naAA) > 2) namsg <- paste0("map_taxa", basetxt, ": can't map groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%), ", ordergroups[3], " (", orderpercent[3], "%)")
    if(sum(naAA) > 3) namsg <- paste0("map_taxa", basetxt, ": can't map groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%), ", sum(naAA) - 2, " others (", sum(orderpercent[-(1:2)]), "%)")
    if(sum(naAA) > 0) print(namsg)
    # Print message about total mapped percent 20200927
    mappedpercent <- formatC(100 - sum(napercent), 1, format = "f")
    print(paste0("map_taxa", basetxt, ": mapping rate to ", refdb, " taxonomy is ", mappedpercent, "%"))
  }
  # Set attributes to indicate unmapped groups 20211007
  attr(iAA, "unmapped_groups") <- nagroups
  attr(iAA, "unmapped_percent") <- napercent
  # Return result
  iAA

}
