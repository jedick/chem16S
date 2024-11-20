# Commands modified from https://benjjneb.github.io/dada2/tutorial.html

library(dada2); packageVersion("dada2")
## [1] ‘1.32.0’

# Files downloaded from https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
##                               reads.in reads.out
## F3D0_S188_L001_R1_001.fastq       7793      7113
## F3D141_S207_L001_R1_001.fastq     5958      5463
## F3D142_S208_L001_R1_001.fastq     3183      2914
## F3D143_S209_L001_R1_001.fastq     3178      2941
## F3D144_S210_L001_R1_001.fastq     4827      4312
## F3D145_S211_L001_R1_001.fastq     7377      6741

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plotErrors(errF, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

# Construct sequence table:w
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
## [1]  20 293

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 251 252 253 254 255 
##   1  88 196   6   2 

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
## [1]  20 232
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.9640374

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#head(track)

# Assign taxonomy
#taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "GTDB_bac120_arc53_ssu_r220_genus.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
##      Kingdom    Phylum         Class         Order           Family           Genus
## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae" NA
## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae" NA
## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae" "UBA7173"
## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae" "UBA7173"
## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae" "Bacteroides"
## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae" "UBA7173"


# Bonus: Handoff to phyloseq
library(phyloseq); packageVersion("phyloseq")
## [1] ‘1.48.0’
library(Biostrings); packageVersion("Biostrings")
## [1] ‘2.72.1’
#library(ggplot2); packageVersion("ggplot2")
#theme_set(theme_bw())

# Construct a simple sample data.frame from the information encoded in the filenames
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

# Store the DNA sequences of our ASVs in the refseq slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

mouse.GTDB_220 <- ps
save(mouse.GTDB_220, file = "mouse.GTDB_220.rda")


# Visualize alpha-diversity:
#plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Ordinate:
# Transform data to proportions as appropriate for Bray-Curtis distances
#ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
