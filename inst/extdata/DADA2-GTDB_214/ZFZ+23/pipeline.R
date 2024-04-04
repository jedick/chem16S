# DADA2/ZFZ+23/pipeline.R
# Pipeline for analyzing Zhang et al. (2023) data with DADA2
# 20230628 Adapted by Jeffrey Dick from the DADA2 Pipeline Tutorial

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

# -------- Sequence processing
# This is adapted from the DADA2 Pipeline Tutorial (version 1.16)
# (https://benjjneb.github.io/dada2/tutorial.html)

###
### Getting ready
###

library(dada2); packageVersion("dada2")
## [1] ‘1.28.0’

path = "."

# Forward fastq filenames have format: SAMPLENAME_1.fastq
# (Reverse reads are only available for two runs in this dataset, so we don't use them)
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###
### Inspect read quality profiles
###

plotQualityProfile(fnFs[1:2])
savePlot("plotQualityProfile.png")
## Plot shows high quality at nearly all positions

###
### Filter and trim
###

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, truncLen = 244,
  maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE)
head(out)
##                     reads.in reads.out
## SRR20752857_1.fastq   138373    129178
## SRR20752858_1.fastq   131800    125352
## SRR20752859_1.fastq   137258    126818
## SRR20752860_1.fastq   149274    141049
## SRR20752861_1.fastq   146572    133904
## SRR20752862_1.fastq   148675    135585

###
### Learn the Error Rates
###

# Timing on 2018 8-core Intel laptop
system.time(errF <- learnErrors(filtFs, multithread = TRUE))
## 127464868 total bases in 522397 reads from 4 samples will be used for learning the error rates.
##      user    system   elapsed
## 10404.745    40.656  1882.410

#saveRDS(errF, "errF.rds")

# Visualize the estimated error rates:
plotErrors(errF, nominalQ = TRUE)
savePlot("plotErrors.png")

###
### Sample Inference
###

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
## Sample 1 - 129178 reads in 18642 unique sequences.
## Sample 2 - 125352 reads in 16697 unique sequences.
## Sample 3 - 126818 reads in 42730 unique sequences.
## Sample 4 - 141049 reads in 40284 unique sequences.
## Sample 5 - 133904 reads in 31922 unique sequences.
## Sample 6 - 135585 reads in 41446 unique sequences.
## Sample 7 - 134653 reads in 60405 unique sequences.

# Inspecting the returned dada-class object:
dadaFs[[1]]
## dada-class: object describing DADA2 denoising results
## 586 sequence variants were inferred from 18642 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

###
### Construct sequence table
###

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
## [1]     7 12439
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
##   244
## 12439

#saveRDS(seqtab, "seqtab.rds")

###
### Remove chimeras
###

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
## [1]    7 9466
sum(seqtab.nochim) / sum(seqtab)
## [1] 0.8360862

saveRDS(seqtab.nochim, "seqtab.nochim.rds")

###
### Track reads through the pipeline
###

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
track
##              input filtered denoisedF nonchim
## SRR20752857 138373   129178    127935  111962
## SRR20752858 131800   125352    124254  120784
## SRR20752859 137258   126818    119463  104525
## SRR20752860 149274   141049    136104  123186
## SRR20752861 146572   133904    129233  101568
## SRR20752862 148675   135585    126366   94678
## SRR20752863 144263   134653    122995   84362

###
### Assign taxonomy
###

# GTDB r214 16S rRNA sequences formatted for DADA2 were downloaded from https://zenodo.org/record/10403693
system.time(taxa <- assignTaxonomy(seqtab.nochim, "/home/sequence/DADA2/GTDB_bac120_arc53_ssu_r214_genus.fa.gz", multithread = TRUE, tryRC = TRUE))
##      user    system   elapsed
## 10143.656     6.935  1312.774

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
##      Kingdom    Phylum             Class             Order               Family            Genus
## [1,] "Bacteria" "Aquificota"       "Aquificae"       "Aquificales"       "Aquificaceae"    "UBA11096"
## [2,] "Bacteria" "Aquificota"       "Aquificae"       "Aquificales"       "Aquificaceae"    "UBA11096"
## [3,] "Bacteria" "Campylobacterota" "Campylobacteria" "Campylobacterales" "Arcobacteraceae" "Aliarcobacter"
## [4,] "Bacteria" "Aquificota"       "Aquificae"       "Aquificales"       "Aquificaceae"    "UBA11096"
## [5,] "Bacteria" "Aquificota"       "Aquificae"       "Aquificales"       "Aquificaceae"    "UBA11096"
## [6,] "Bacteria" "Campylobacterota" "Campylobacteria" "Campylobacterales" "Arcobacteraceae" "Aliarcobacter"

saveRDS(taxa, "taxa.rds")

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

# Read data frame with sample data constructed from SRA entries and Table S5 of Zhang et al. (2023)
samdf <- read.csv("sample_data.csv", row.names = 1)
# Make sure sample names of sequences and sample data are the same
stopifnot(all(rownames(seqtab.nochim) == rownames(samdf)))
# We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
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

saveRDS(ps, "ps_ZFZ+23.rds", compress = "xz")

# Visualize alpha-diversity:
plot_richness(ps, x = "SampleName", measures = c("Shannon", "Simpson"), color = "T")
savePlot("plot_richness.png")

# Ordinate:
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "T", title = "Bray NMDS")
savePlot("plot_ordination.png")

# Bar plot:
top100 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x = "T", fill = "Class")
savePlot("plot_bar.png")

# -------- End of sequence processing pipeline
