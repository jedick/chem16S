# dada2/FEN+22/pipeline.R
# Pipeline for analyzing Fonseca et al. (2022) data with DADA2
# 20230625 jmd

# -------- Prepare FASTQ files
# Read SRA run info file from NCBI
sra <- read.csv("SraRunInfo.csv")
# Create FASTQ files for each run
for(run in sra$Run) {
  cmd <- paste("fastq-dump --split-files --skip-technical --clip", run)
  print(cmd)
  system(cmd)
}
# --------

# -------- Notes for processing 454 data
# (https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data)

###
### Can I use dada2 with my 454 or Ion Torrent data?
###

## Yes. We recommend the following parameters for denoising pyrosequencing data (like IT and 454):
## dada(..., HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

## We recommend additional (relative to the tutorial) filtering of 454 sequences by maximum length:
## filterAndTrim(..., maxLen=XXX) # XXX depends on the chemistry
# --------

# -------- Sequence processing
# This is adapted from the DADA2 pipeline tutorial (version 1.16)
# (https://benjjneb.github.io/dada2/tutorial.html)

###
### Getting ready
###

library(dada2); packageVersion("dada2")
## [1] ‘1.28.0’

path = "."

# Forward fastq filenames have format: SAMPLENAME_2.fastq
# (This dataset doesn't have reverse reads)
fnFs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###
### Inspect read quality profiles
###

plotQualityProfile(fnFs[1:2])
savePlot("plotQualityProfile.png")
## Plot shows quality-score dropoff at ~500 nt

###
### Filter and trim
###

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

# Notes:
# maxEE = 5 is increased from setting in tutorial (2) because 454 has longer sequences than Illumina
# maxLen = 550 is to filter sequences by maximum length
# Truncation to 500 nt is performed *after* maxLen filtering (see ?filterAndTrim)
out <- filterAndTrim(fnFs, filtFs, maxLen = 550, truncLen = 500,
  maxN = 0, maxEE = 5, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE)
head(out)
##                    reads.in reads.out
## SRR1346068_2.fastq     9102      1539
## SRR1346069_2.fastq     9756      1641
## SRR1346070_2.fastq     6876       727
## SRR1346071_2.fastq     6573      1018
## SRR1346072_2.fastq    10375      1855
## SRR1346073_2.fastq    10505      1769

###
### Learn the Error Rates
###

# Timing on 2018 Intel notebook
system.time(errF <- learnErrors(filtFs, multithread = TRUE))
##     user   system  elapsed
## 4406.295    6.821  694.335

#saveRDS(errF, "errF.rds")

# Visualize the estimated error rates:
plotErrors(errF, nominalQ = TRUE)
savePlot("plotErrors.png")

###
### Sample Inference
###

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
## Sample 1 - 1539 reads in 1535 unique sequences.
## Sample 2 - 1641 reads in 1641 unique sequences.
## Sample 3 - 727 reads in 727 unique sequences.
## Sample 4 - 1018 reads in 1018 unique sequences.
## Sample 5 - 1855 reads in 1855 unique sequences.
## Sample 6 - 1769 reads in 1767 unique sequences.
## Sample 7 - 1951 reads in 1951 unique sequences.
## Sample 8 - 3236 reads in 3190 unique sequences.
## Sample 9 - 2059 reads in 2027 unique sequences.
## Sample 10 - 3011 reads in 2959 unique sequences.
## Sample 11 - 4171 reads in 4112 unique sequences.
## Sample 12 - 2258 reads in 2231 unique sequences.
## Sample 13 - 2444 reads in 2416 unique sequences.
## Sample 14 - 2074 reads in 2044 unique sequences.
## Sample 15 - 18050 reads in 17955 unique sequences.
## Sample 16 - 26949 reads in 26755 unique sequences.
## Sample 17 - 6170 reads in 6106 unique sequences.
## Sample 18 - 7681 reads in 7626 unique sequences.
## Sample 19 - 21516 reads in 21408 unique sequences.
## Sample 20 - 18801 reads in 18699 unique sequences.
## Sample 21 - 9287 reads in 9206 unique sequences.
## Sample 22 - 13607 reads in 12664 unique sequences.
## Sample 23 - 12567 reads in 12058 unique sequences.
## Sample 24 - 14769 reads in 13744 unique sequences.
## Sample 25 - 12823 reads in 12080 unique sequences.
## Sample 26 - 9773 reads in 9271 unique sequences.
## Sample 27 - 12838 reads in 11773 unique sequences.
## Sample 28 - 14316 reads in 10066 unique sequences.

# Inspecting the returned dada-class object:
dadaFs[[1]]
## dada-class: object describing DADA2 denoising results
## 4 sequence variants were inferred from 1535 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

###
### Construct sequence table
###

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
## [1]   28 2546
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
##  500
## 2546

#saveRDS(seqtab, "seqtab.rds")

###
### Remove chimeras
###

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
## [1]   28 2265
sum(seqtab.nochim) / sum(seqtab)
## [1] 0.8556967

#saveRDS(seqtab.nochim, "seqtab.nochim.rds")

###
### Track reads through the pipeline
###

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track
##            input filtered denoisedF nonchim
## SRR1346068  9102     1539        39      39
## SRR1346069  9756     1641        13      13
## SRR1346070  6876      727        65      65
## SRR1346071  6573     1018         1       1
## SRR1346072 10375     1855         7       7
## SRR1346073 10505     1769         4       4
## SRR1346074  9890     1951         6       6
## SRR1346075 22382     3236       626     626
## SRR1346076 13661     2059       360     360
## SRR1346077 19838     3011       698     694
## SRR1346078 22505     4171       563     563
## SRR1346079 17407     2258       139     139
## SRR1346080 17975     2444        96      96
## SRR1346081 14498     2074       187     187
## SRR1346082 29397    18050      2324    2204
## SRR1346083 40208    26949      3917    3146
## SRR1346084 10887     6170      1807    1760
## SRR1346085 13557     7681      1782    1742
## SRR1346086 34193    21516      3753    2776
## SRR1346087 32586    18801      3384    3120
## SRR1346088 13847     9287      1353    1180
## SRR1346089 16638    13607      6450    5539
## SRR1346090 17963    12567      4020    3734
## SRR1346091 19135    14769      5764    4912
## SRR1346092 18719    12823      5429    5078
## SRR1346093 16122     9773      3355    2973
## SRR1346094 19837    12838      6243    5219
## SRR1346095 15701    14316     10739    7832

# NOTE: There are low numbers of sequences in Batch A (first 14 samples)
# Filter seqtab.nochim to keep only Batch B (last 14 samples)
seqtab.nochim.B <- seqtab.nochim[15:28, ]
# Some ASVs have zero counts in the table ...
sum(colSums(seqtab.nochim.B) == 0)
## [1] 170
# ... so let's remove them
seqtab.nochim.B <- seqtab.nochim.B[, colSums(seqtab.nochim.B) > 0]
dim(seqtab.nochim.B)
## [1]   14 2095

#saveRDS(seqtab.nochim.B, "seqtab.nochim.B.rds")

###
### Assign taxonomy
###

# The GTDB r207 training set for DADA2 is available at https://zenodo.org/record/6655692
system.time(taxa <- assignTaxonomy(seqtab.nochim.B, "/home/sequence/dada2/GTDB_bac120_arc53_ssu_r207_Genus.fa.gz", multithread = TRUE))
##     user   system  elapsed
## 3873.874    6.198  541.574

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
##      Kingdom    Phylum           Class                 Order              Family              Genus
## [1,] "Bacteria" "Bacteroidota"   "Bacteroidia"         "Flavobacteriales" "Flavobacteriaceae" NA
## [2,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" NA                 NA                  NA
## [3,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Beggiatoales"     "Beggiatoaceae"     "Marithrix"
## [4,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Beggiatoales"     "Beggiatoaceae"     "Marithrix"
## [5,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Beggiatoales"     "Beggiatoaceae"     "Marithrix"
## [6,] "Bacteria" "Bacteroidota"   "Bacteroidia"         "Flavobacteriales" "Flavobacteriaceae" "Flavicella"

#saveRDS(taxa, "taxa.rds")

###
### Bonus: Handoff to phyloseq
###

library(phyloseq); packageVersion("phyloseq")
## [1] ‘1.44.0’
library(Biostrings); packageVersion("Biostrings")
## [1] ‘2.68.1’
library(ggplot2); packageVersion("ggplot2")
## [1] ‘3.4.2’
theme_set(theme_bw())

# Read data frame with sample data constructed from BioSample and SRA entries
samdf <- read.csv("sample_data.csv", row.names = 1)
# Make sure sample names of sequences and sample data are the same
stopifnot(all(rownames(seqtab.nochim.B) == rownames(samdf)))
# We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim.B, taxa_are_rows = FALSE),
  sample_data(samdf), tax_table(taxa))

# Use short taxa names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2095 taxa and 14 samples ]
## sample_data() Sample Data:       [ 14 samples by 19 sample variables ]
## tax_table()   Taxonomy Table:    [ 2095 taxa by 6 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 2095 reference sequences ]

saveRDS(ps, "ps_FEN+22.rds", compress = "xz")

# Visualize alpha-diversity:
plot_richness(ps, x = "Location", measures = c("Shannon", "Simpson"), color = "Sediment_redox")
savePlot("plot_richness.png")

# Ordinate:
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "Location", title = "Bray NMDS")
savePlot("plot_ordination.png")

# One sample is clearly different from the others.
# It is from 50 m depth at Valparaiso sampled on 2012-05-12.
# Note that this sample is present in SRA but was not analyzed by Fonseca et al. (2012)
#   (their Table 1 has 13 of the 14 available samples).
# Let's remove this sample and try the ordination again.
ps <- prune_samples(sample_names(ps) != "SRR1346095", ps)
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
p2 <- plot_ordination(ps.prop, ord.nmds.bray, color = "Location", title = "Bray NMDS")
p2 + geom_polygon(aes(fill = Location)) + geom_point(size = 4)
savePlot("plot_ordination2.png")

# Bar plot:
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x = "Location", fill = "Class")
savePlot("plot_bar.png")


# -------- End of sequence processing pipeline
