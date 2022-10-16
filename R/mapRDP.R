# chem16S/mapRDP.R
# Map RDP to RefSeq taxonomy 20200912
# Moved to chem16S 20220505
# Add refdb argument 20221016

mapRDP <- function(RDP = NULL, refdb = "RefSeq", quiet = TRUE) {

  # Make group names by combining rank and name
  RDPgroups <- paste(RDP$rank, RDP$name, sep = "_")
  # Calculate group abundances for displaying messages
  groupcounts <- as.numeric(rowSums(RDP[, -(1:4), drop = FALSE]))
  # Basename of file to use in messages 20220505
  basefile <- basename(attr(RDP, "file"))

  # Manual mappings for RefSeq (NCBI) taxonomy but not for GTDB 20221016
  if(refdb == "RefSeq") {

    NCBIgroups <- vapply(RDPgroups, switch, "",
      # 20200920 Lots of Escherichia in urine [WZZ+18]
      "genus_Escherichia/Shigella" = "genus_Escherichia",
      "phylum_Cyanobacteria/Chloroplast" = "phylum_Cyanobacteria",
      # 20200929 Unclassified Cyanobacteria are just Cyanobacteria
      "class_Cyanobacteria" = "phylum_Cyanobacteria",
      "genus_Spartobacteria_genera_incertae_sedis" = "class_Spartobacteria",
      # 20210502 Processing Guerrero Negro
      "class_Planctomycetacia" = "class_Planctomycetia",
      # 20210526 NCBI taxonomy no longer has an Actinobacteria "class"
      "class_Actinobacteria" = "phylum_Actinobacteria",
      "order_Rhizobiales" = "order_Hyphomicrobiales",
      # 20210530 Acidobacteria
      "genus_Gp1" = "genus_Acidobacterium",
      "genus_Gp6" = "genus_Luteitalea",
      # 20210530 Cyanobacteria
      "genus_GpI" = "genus_Nostoc",
      "genus_GpIIa" = "genus_Synechococcus",
      "genus_GpVI" = "genus_Pseudanabaena",
      "family_Family II" = "family_Synechococcaceae",
      # 20210609 Verrucomicrobia
      "genus_Subdivision3_genera_incertae_sedis" = "family_Verrucomicrobia subdivision 3",
      # 20211215 Clostridia
      # Clostridiales is a synonym for Eubacteriales (https://lpsn.dsmz.de/order/eubacteriales)
      "order_Clostridiales" = "order_Eubacteriales",
      # Ruminococcaceae is a synonym for Oscillospiraceae (https://lpsn.dsmz.de/family/oscillospiraceae)
      "family_Ruminococcaceae" = "family_Oscillospiraceae",

      ## NOT USED

      # 20200929 Yellowstone [BGPF13]
      # Not used because there's not much information about this group
      #"genus_Armatimonadetes_gp7" = "phylum_Armatimonadetes",

      # 20210530 Marcellus Shale [CHM+14]
      # https://lpsn.dsmz.de/family/arcobacteraceae
      # Not used because this causes a large low-ZC deviation in Blue Hole 100m sample
      #"family_Arcobacteraceae" = "family_Campylobacteraceae",

    NA_character_)
    iswitch <- !is.na(NCBIgroups)
    if(any(iswitch)) {
      # Make the switch
      RDPgroups[iswitch] <- NCBIgroups[iswitch]
      if(!quiet) {
        # Print message(s) about switched names and abundance
        from <- names(NCBIgroups)[iswitch]
        to <- NCBIgroups[iswitch]
        switchcounts <- groupcounts[iswitch]
        switchpercent <- round(switchcounts / sum(groupcounts) * 100, 1)
        # Only print message for mappings of groups at least 0.1% abundant 20200927
        if(any(switchpercent >= 0.1)) {
          print(paste0("mapRDP [", basefile, "]: using the following RDP --> RefSeq (NCBI) mapping(s):"))
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
  AApath <- file.path("extdata", refdb, "taxon_AA.csv.xz")
  AAfile <- system.file(AApath, package = "chem16S")
  taxon_AA <- read.csv(AAfile, as.is = TRUE)
  AAgroups <- paste(taxon_AA$protein, taxon_AA$organism, sep = "_")
  # Do the mapping!
  iAA <- match(RDPgroups, AAgroups)
  # Get percentages of unmapped groups
  naAA <- is.na(iAA)
  nagroups <- RDPgroups[naAA]
  nacounts <- groupcounts[naAA]
  napercent <- nacounts / sum(groupcounts) * 100
  if(!quiet) {
    # Print summary of missing groups
    naorder <- order(napercent, decreasing = TRUE)
    ordergroups <- nagroups[naorder]
    orderpercent <- round(napercent[naorder], 2)
    if(sum(naAA) > 0) namsg <- paste0("mapRDP [", basefile, "]: can't map RDP group ", ordergroups[1], " (", orderpercent[1], "%)")
    if(sum(naAA) > 1) namsg <- paste0("mapRDP [", basefile, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%)")
    if(sum(naAA) > 2) namsg <- paste0("mapRDP [", basefile, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%), ", ordergroups[3], " (", orderpercent[3], "%)")
    if(sum(naAA) > 3) namsg <- paste0("mapRDP [", basefile, "]: can't map RDP groups ", ordergroups[1], " (", orderpercent[1], "%), ",
                                      ordergroups[2], " (", orderpercent[2], "%), ", sum(naAA) - 2, " others (", sum(orderpercent[-(1:2)]), "%)")
    if(sum(naAA) > 0) print(namsg)
    # Print message about total mapped percent 20200927
    mappedpercent <- formatC(100 - sum(napercent), 1, format = "f")
    print(paste0("mapRDP [", basefile, "]: mapped ", mappedpercent, "% of RDP classifications to ", refdb, " taxonomy"))
  }
  # Set attributes to indicate unmapped groups 20211007
  attr(iAA, "unmapped_groups") <- nagroups
  attr(iAA, "unmapped_percent") <- napercent
  # Return result
  iAA

}

