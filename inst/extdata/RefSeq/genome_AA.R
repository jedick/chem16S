# chem16S/RefSeq/genome_AA.R
# Calculate the overall amino acid composition of proteins for each taxid in RefSeq
# 20100704 First version
# 20130922 Deal with WP multispecies accessions (RefSeq 61)
# 20190830 Update for RefSeq 95:
# - Start with blank composition table for all taxa
# - Read files in parallel and add sequences serially to composition tables
# - After parallel reading, combine composition tables

# This script depends on system "join" command
# canprot is needed for read_fasta
library(canprot)
# Set maxcores to 1 here -- we will parallelize across files
thermo("opt$maxcores" = 1)
library(parallel)

## Define global variables

# Generate a blank data frame of amino acid compositions for all taxa 20190830
mk.aadat <- function() {
  # get all of the taxids in the release catalog (microbial proteins) and multispecies accessions
  catalog_unique_taxids <- read.table("unique_taxids.txt")$V1
  WP_unique_taxids <- read.table("WP_unique_taxids.txt")$V1
  unique_taxids <- sort(unique(c(catalog_unique_taxids, WP_unique_taxids)))
  # initialize one row of data frame of amino acid composition
  aa1 <- thermo()$protein[1, ]
  aa1[, 1:4] <- ""
  aa1[, 5:25] <- 0
  # make the number of rows equal to the number of taxids
  aa <- aa1[rep(1, length(unique_taxids)), ]
  # add taxids to "organism" column
  aa$protein <- "refseq"
  aa$organism <- unique_taxids
  aa
}
aadat <- mk.aadat()

## Define functions to read files

# Read seqeunces from a single file into a blank aadat data frame 20190830
read.file <- function(file) {
  # Where to save the results
  outfile <- file.path("csv", gsub(".gz", ".csv", basename(file)))
  # Skip existing files (for restarting after an interrupted run) 20200729
  if(file.exists(outfile)) return()

  # Read the fasta file and find the header lines
  lines <- readLines(file)
  ihead <- grep("^>", lines)
  # Get the amino acid compositions from the file
  aa <- suppressMessages(read_fasta(file = "", lines = lines, ihead = ihead))
  # Sort the data by accession (needed for join command)
  accord <- order(aa$protein)
  ihead <- ihead[accord]
  aa <- aa[accord, ]

  # Identify any multispecies (WP) accessions
  # Use unique tempfiles to prevent clashes between parallel processes
  # (but not R's tempfile(), which seems prone to connection errors in parallel runs)
  acc.txt <- paste0(basename(file), ".acc.txt")
  WP.taxid.match <- paste0(basename(file), ".WP.taxid.match")
  write.table(aa$protein, acc.txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Use tab separator here: https://stackoverflow.com/questions/1722353/unix-join-separator-char
  system(paste("join -j 1 -t $'\t'", acc.txt, "release206.MultispeciesAutonomousProtein2taxname >", WP.taxid.match))
  # Turn off quoting because Synechococcus sp. JA-2-3B'a(2-13) messes things up 20190831
  # Include check for empty file (no multispecies accessions) 20190901
  if(file.size(WP.taxid.match) > 0) myWP <- read.table(WP.taxid.match, sep = "\t", stringsAsFactors = FALSE, quote = "")
  else myWP <- data.frame(V1 = character(), V2 = numeric())
  message(paste0(file, " has ", nrow(aa), " sequences (", length(unique(myWP$V1)),
    " multispecies accessions for ", length(unique(myWP$V2)), " taxids)"))
  file.remove(WP.taxid.match)

  # Get the organism names
  organism <- gsub("\\]$", "", gsub(".*\ \\[", "", lines[ihead]))
  # Find the taxids that are associated with these accessions
  acc.taxid.match <- paste0(basename(file), ".acc.taxid.match")
  system(paste("join -j 1", acc.txt, "accession.taxid.txt >", acc.taxid.match))
  acc.taxid <- read.table(acc.taxid.match)
  file.remove(acc.txt)
  file.remove(acc.taxid.match)

  # Create a temporary data frame to hold the aggregated amino acid counts
  aatmp <- aa[, c(2, 5:25)]
  # Put the taxid in the "organism" column
  aatmp$organism <- acc.taxid$V2
  # Aggregate (sum) the amino acid composition on taxid
  # https://stackoverflow.com/questions/34523679/aggregate-multiple-columns-at-once
  aaagg <- aggregate(.~organism, aatmp, sum)
  # Get organism names for each taxid
  iacc <- match(aaagg$organism, acc.taxid$V2)
  orgnames <- organism[iacc]
  # Get all matching taxids
  iaa <- match(aaagg$organism, aadat$organism)
  # Add amino acid composition and organism name for all matching taxids
  aadat[iaa, 5:25] <- aadat[iaa, 5:25] + aaagg[, 2:22]
  aadat$ref[iaa] <- orgnames

  # Check for multispecies accessions
  if(nrow(myWP) > 0) {
    isWP <- aa$protein %in% myWP$V1
    # Loop over the sequences to add amino acid compositions
    for(i in which(isWP)) {
      iWP <- myWP$V1 == aa$protein[i]
      taxid <- myWP$V2[iWP]
      iaa <- match(taxid, aadat$organism)
      # Repeat the amino acid composition by the number of matching WP taxids
      aa.WP <- aa[rep(i, sum(iWP)), 5:25]
      # Add the amino acid composition and organism names for the matching taxids
      aadat[iaa, 5:25] <- aadat[iaa, 5:25] + aa.WP
      aadat$ref[iaa] <- myWP$V3[iWP]
    }
  }

  # Save CSV file
  write.csv(aadat, outfile, row.names = FALSE, quote = 3)
  return()
}

read.allfiles <- function() {
  # List all files in "protein" directory
  files <- file.path("protein", dir("protein"))
  # Exclude already processed files 20210604
  csvfiles <- file.path("csv", dir("csv"))
  done <- gsub("csv$", "gz", gsub("^csv", "protein", csvfiles))
  files <- files[!files %in% done]

  ## Parallel computations, adapted from ?clusterApply
  # Use option cl.cores to choose an appropriate cluster size.
  # Use type = "FORK" to share address space/environment: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
  # -- No, don't use it, because of file connection problems
  # Use outfile = "" to show messages: https://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
  cl <- makeCluster(getOption("cl.cores", 8), outfile = "")
  # Set up R environment (would not be needed with FORK cluster)
  clusterExport(cl, c("read.file", "aadat"))
  clusterEvalQ(cl, suppressMessages(library(canprot)))

  ## For debugging: https://stackoverflow.com/questions/16895848/results-of-workers-not-returned-properly-snow-debug
  #workerfun <- function(i) {
  #tryCatch({
  #  read.file(i)
  #},
  #error=function(e) {
  #  print(e)
  #  stop(e)
  #})
  #}
  #aalist <- parLapply(cl, files, workerfun)

  # Save a CSV file of amino acid compositions (per taxid) for each fasta file
  parLapply(cl, files, read.file)
  stopCluster(cl)
  #invisible(aalist)
}

genome_AA <- function() {
  # List all files in "csv" directory
  files <- file.path("csv", dir("csv"))
  # Read the CSV files into a single list
  aalist <- lapply(files, read.csv, as.is = TRUE)
  # Extract the numeric parts of the data frames to be added
  aamat <- lapply(aalist, function(x) x[, 5:25])
  # Add together the amino acid compositions from all files
  sumaa <- Reduce("+", aamat)
  # Combine the organism names from all the files
  aaref <- sapply(aalist, function(x) x$ref)
  # Prepend an all-"" column first so that "" always appears first in unique()
  aaref <- cbind(aadat$ref, aaref)
  # Get unique names, separated by ";"
  aaref <- sapply(apply(aaref, 1, unique), function(x) paste(x, collapse = ";"))
  # Remove the leading ;
  aaref <- gsub("^;", "", aaref)
  # Make a combined data frame
  aaall <- cbind(aadat[, 1:4], sumaa)
  aaall$ref <- aaref
  # Save the result
  write.csv(aaall, "genome_AA.csv", row.names = FALSE, quote = 3)
  invisible(aaall)
}
