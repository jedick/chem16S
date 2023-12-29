# Generate a table with names of
# species,genus,family,order,class,phylum,superkingdom
# for species-level taxa in RefSeq database

# Uses functions in CHNOSZ to process taxonomy files
require(CHNOSZ)

# Change this to the location where names.dmp and nodes.dmp are located
taxdir <- "./taxdump"

# Get the taxids from genome_AA.csv
pr <- read.csv("genome_AA.csv")
taxid <- pr$organism

# Read in the names and nodes
if(!exists("taxnames")) {
  cat("reading names...\n")
  taxnames <- getnames(taxdir)
}
if(!exists("taxnodes")) {
  cat("reading nodes...\n")
  taxnodes <- getnodes(taxdir)
}

# What ranks we want to get
ranks <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")

# Start with an empty list
out <- rep(list(character()), length(ranks))
names(out) <- c(ranks)

# Loop over taxids
ii <- seq_along(taxid)
for(i in ii) {
  # Test if data are available for this taxid
  if(!taxid[i] %in% taxnames$id) {
    for(j in 1:length(ranks)) out[[j]] <- c(out[[j]],NA)
    next
  }
  # Get taxids of all parents
  pids <- allparents(taxid[i], taxdir, nodes = taxnodes)
  # Get ranks of all parents
  pranks <- getrank(pids, taxdir, nodes = taxnodes)
  # Find which parents are in the required ranks
  ip <- match(ranks, pranks)
  # Get names of these parents
  pnames <- sciname(pids[ip], taxdir, names = taxnames)
  # Add results to output list
  for(j in 1:length(ranks)) out[[j]] <- c(out[[j]], pnames[j])
  # Report progress
  if(i %% 50 == 0) cat(paste(i, ".. "))
}
# Rinish progress report
cat("done!\n")

# Write results to a file
out <- as.data.frame(out)
out <- cbind(data.frame(taxid = taxid[ii], out))
write.csv(out, "taxonomy.csv", row.names = FALSE, quote = FALSE)
