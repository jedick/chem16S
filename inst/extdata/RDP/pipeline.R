# chem16S/RDP/pipeline.R
# Merge and filter 16S sequences from SRA or MG-RAST, subsample, find chimeras, and run RDP Classifier

# 20200909 Initial version by Jeffrey M. Dick
# 20210922 Add settings for orp16S datasets
#          Add --skip-technical --clip options to fastq-dump
#          Don't use fastq_filter --stripleft for 454 data

## SYSTEM REQUIREMENTS
# (version numbers for information only)

# R 4.0.0              https://www.r-project.org/
# fastq-dump 2.9.0     https://ncbi.github.io/sra-tools/
# vsearch 2.15.0       https://github.com/torognes/vsearch
# seqtk 1.3-r115       https://github.com/lh3/seqtk
# RDP Classifier 2.13  https://sourceforge.net/projects/rdp-classifier/
# Java openjdk version 1.8.0_252 (for RDP); GNU utils (cat, awk, grep)

## USAGE

# First, make sure the working directory is clean (especially no *.fa files)
# To filter one file (generates .fa files in workdir for subsampling)
# filter("SRR2059382")

# To filter multiple files:
# RUNID <- c("SRR2059381", "SRR2059380", "SRR2059379")
# lapply(RUNID, filter)

# To subsample and remove chimeras from the filtered files:
# subsample()
# findchimeras()

# To run RDP Classifier after chimera removal:
# classify(RUNID)
# To combine output of RDP Classifier into one file
# -- result is saved in RDP/<study>.tab in the current directory
# mkRDP(study, RUNID)

## STUDY SETTINGS

# Change the following line to setup the pipeline for one study
study <- "SMS+12"
# Settings for all studies are stored here
file <- tempfile()
# Write spaces here (but don't save them) to make this easier to read
writeLines(con = file, text = gsub(" ", "", c(
  # For 454 experiments, set forwardonly to NA
  "study, forwardonly, trunclen",
  # For HCW+13 (Guerrero Negro), don't use filter(); start with findchimeras()
  #"HCW+13, NA, NA",     # GenBank JN427016-JN539989  10.1038/ismej.2012.79
  "BGPF13, NA, NA",      # SRA PRJNA207095            10.3389/fmicb.2013.00330
  "HLA+16, NA, NA"       # SRA PRJEB1245              10.3389/fmicb.2016.01883
  # For SMS+12 (Bison Pool), start with classify()
  #"SMS+12, NA, NA"      # MG-RAST mgp100292          10.1371/journal.pone.0038108
)))

# This reads and applies the settings
settings <- read.csv(file, as.is = TRUE)
istudy <- match(study, settings$study)
# Flag for 454 studies
is454 <- is.na(settings$forwardonly[istudy])
# Use forward reads only (for low-quality or missing reverse reads)
forwardonly <- settings$forwardonly[istudy]
# Value of fastq_trunclen
trunclen <- settings$trunclen[istudy]

## FILE LOCATIONS

# Working directory
workdir <- "~/tmp/chem16S"
# Output directory for FASTA files
FASTAdir <- file.path("/home/chem16S", study, "fasta")
# Reference database for chimera detection
refdb <- "/home/chem16S/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
#refdb <- "/home/chem16S/chimera/gold.fa"
# RDP jar file
RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
# Output directory for RDP Classifier
RDPdir <- file.path("/home/chem16S", study, "RDP")
# RDP Classifier version (java -jar /opt/rdp_classifier_2.13/dist/classifier.jar version)
#Gene:16srrna    Trainset No:18  Taxonomy Version:RDP 16S rRNA training setNo 18 07/2020
#RDP Classifier Version:RDP Naive Bayesian rRNA Classifier Version 2.11, September 2015

filter <- function(RUNID) {
  message("============================================")
  message(paste0("Merging and filtering ", RUNID, " [", study, "]"))
  message("============================================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))
  # The output file from this function is a FASTA file with .fa suffix
  outfile <- paste0(RUNID, ".fa")

  if(study %in% c("BYB+17", "GRG+20", "MKK+11", "FLA+12")) {
    # For BYB+17, FASTQ files are downloaded from SRA cloud 20210917
    # For the others, files are downloaded from MG-RAST 20220122
    fqdump <- FALSE
  } else {
    # Generate input FASTQ files with fastq-dump
    fqdump <- TRUE
    cmd <- paste("fastq-dump --split-files --skip-technical --clip", RUNID)
    if(study == "SCH+16") cmd <- paste("fastq-dump --split-files --skip-technical --clip -X 100000", RUNID)
    print(cmd)
    system(cmd)
  }

  if(study %in% c("GRG+20", "MKK+11", "FLA+12")) {
    # These are files from MG-RAST (.fasta suffix)
    # Prefix MG-RAST ID to header so we can extract the reads in the subsample() and findchimeras() steps 20211004
    fastafile <- paste0(RUNID, ".fasta")
    lines <- readLines(fastafile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(">", RUNID, ";", gsub("^>", "", lines[ihead]))
    writeLines(lines, outfile)
  } else if(is454) {
    # For 454 studies, fastq-dump --skip-technical puts biological sequences into file with highest-numbered suffix 
    # Use _3.fastq, _2.fastq or _1.fastq if they are present 20200921
    infile <- paste0(RUNID, "_3.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_2.fastq")
    if(!file.exists(infile)) infile <- paste0(RUNID, "_1.fastq")
    # This is needed for files not generated by fastq-dump (BYB+17)
    if(!file.exists(infile)) infile <- paste0(RUNID, ".fastq")
    truncqual <- 15
    fastq_minlen <- 200
    if(study %in% c("SAR+13", "GFE+16")) truncqual <- 11
    if(study %in% c("GFE+16")) fastq_minlen <- 100
    cmd <- paste("vsearch -fastq_filter", infile, "-fastq_minlen", fastq_minlen, "-fastq_maxlen 600 -fastq_truncqual", truncqual, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  } else if(forwardonly) {
    # Use forward reads only
    file.copy(paste0(RUNID, "_1.fastq"), "merged.fastq", overwrite = TRUE)
    # For Human Microbiome Project, reads are in different file 20211216
    if(study == "HMP12") file.copy(paste0(RUNID, "_4.fastq"), "merged.fastq", overwrite = TRUE)
    nseq <- length(readLines("merged.fastq")) / 4
    print(paste0("Using forward reads only (", nseq, " sequences)"))
  } else {
    # Check that we have same number of forward and reverse reads 20200917
    # (issue with some datasets: fastq_dump rejects some reads because READLEN < 1)
    ID1 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", RUNID, "_1.fastq"), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", RUNID, "_2.fastq"), intern = TRUE)), " "), "[", 1)
    # Forward reads that have matching reverse reads
    forward_matches <- ID1[ID1 %in% ID2]
    if(length(forward_matches) < length(ID1)) {
      writeLines(forward_matches, "forward_matches.txt")
      print(cmd <- paste0("seqtk subseq ", RUNID, "_1.fastq forward_matches.txt > new_1.fastq"))
      system(cmd)
      print(cmd <- paste0("mv new_1.fastq ", RUNID, "_1.fastq"))
      system(cmd)
    }
    # Reverse reads that have matching forward reads
    reverse_matches <- ID2[ID2 %in% ID1]
    if(length(reverse_matches) < length(ID2)) {
      writeLines(reverse_matches, "reverse_matches.txt")
      print(cmd <- paste0("seqtk subseq ", RUNID, "_2.fastq reverse_matches.txt > new_2.fastq"))
      system(cmd)
      print(cmd <- paste0("mv new_2.fastq ", RUNID, "_2.fastq"))
      system(cmd)
    }
    # Make sure forward and reverse reads have same IDs
    ID1 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", RUNID, "_1.fastq"), intern = TRUE)), " "), "[", 1)
    ID2 <- sapply(strsplit(gsub("@", "", system(paste0("awk '(NR - 1) % 4 == 0' ", RUNID, "_2.fastq"), intern = TRUE)), " "), "[", 1)
    if(!all(ID1 == ID2)) stop(paste0("sequence IDs in ", RUNID, "_1.fastq and ", RUNID, "_2.fastq don't match!"))

    # Merge paired reads
    print(cmd <- paste0("vsearch -fastq_mergepairs ", RUNID, "_1.fastq -reverse ", RUNID, "_2.fastq -fastqout merged.fastq"))
    system(cmd)
  }
  if(!is454) {
    # Remove short and high-error reads
    # Use fastq_maxee_rate instead of fastq_maxee 20200920
    maxee_rate <- 0.005
    if(study == "YHK+19") maxee_rate <- 0.05
    # Use default -fastq_qmax except for some studies 20210922
    qmax <- 41
    # For Illumina MiSeq (strange quality encoding?) 20210909
    if(study == "ECS+18") qmax <- 93
    # For Ion Torrent 20210928
    if(study == "PCL+18") qmax <- 45
    cmd <- paste("vsearch -fastq_filter merged.fastq -fastq_trunclen", trunclen, "-fastq_qmax", qmax, "-fastq_maxee_rate", maxee_rate, "-fastaout", outfile)
    print(cmd)
    system(cmd)
  }
  if(study %in% c("BYB+17")) {
    # Append run ID to header so we can extract the reads in the subsample() and findchimeras() steps 20210913
    lines <- readLines(outfile)
    ihead <- grep("^>", lines)
    lines[ihead] <- paste0(lines[ihead], ";", RUNID)
    writeLines(lines, outfile)
  }

  # Clean up
  if(fqdump) file.remove(Sys.glob("*.fastq"))
  return()
}

subsample <- function(do.singletons = TRUE) {
  message("=======================")
  message(paste0("Subsampling [", study, "]"))
  message("=======================")

  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))

  if(do.singletons) {
    # Make combined FASTA file of filtered sequences 20200914
    print(cmd <- "cat *.fa > filtered.fasta")
    system(cmd)
    # Remove singletons 20200912
    # 1. Dereplicate sequences, keeping those that appear exactly once
    print(cmd <- "vsearch -derep_fulllength filtered.fasta -maxuniquesize 1 -output singletons.fasta")
    system(cmd)
    # 2. Get all sequence IDs
    print(cmd <- 'grep ">" filtered.fasta | sed -e s/\\ .*//g | sed -e s/\\>//g > allIDs.txt')
    system(cmd)
    allIDs <- readLines("allIDs.txt")
    # 3. Get singleton sequence IDs
    print(cmd <- 'grep ">" singletons.fasta | sed -e s/\\ .*//g | sed -e s/\\>//g > singletonIDs.txt')
    system(cmd)
    singletonIDs <- readLines("singletonIDs.txt")
    # 4. Write IDs of sequences that are not singletons
    writeLines(setdiff(allIDs, singletonIDs), "not_singletons.txt")
    # 5. Extract sequences that are not singletons
    print(cmd <- "seqtk subseq filtered.fasta not_singletons.txt > not_singletons.fasta")
    system(cmd)
    singletxt <- paste0("Removed ", length(singletonIDs), " singletons from ", length(allIDs), " sequences (", round(length(singletonIDs) / length(allIDs) * 100, 1), "%)")
    print(singletxt)
  } else {
    # Don't remove singletons, but combine filtered sequences for subsampling 20210821
    print(cmd <- "cat *.fa > not_singletons.fasta")
    system(cmd)
  }

  # Loop over samples
  allRUNID <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisRUNID in allRUNID) {
    # Extract sequences for this sample
    thisfile <- paste0(thisRUNID, ".fasta")
    # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
    # Don't use "awk '/>" because sample names for SBE+17 are at end of header of original FASTQ files, not beginning 20210501
    print(cmd <- paste0("awk '/", thisRUNID, "/{p++;print;next} /^>/{p=0} p' not_singletons.fasta > ", thisfile))
    system(cmd)

    # Subsample 10000 sequences - put output in .fa file in FASTAdir
    print(cmd <- paste('grep "^>"', thisfile))
    thisfile.headers <- system(cmd, intern = TRUE)
    nout <- nseq <- length(thisfile.headers)
    outfile <- file.path(FASTAdir, paste0(thisRUNID, ".fa"))
    if(nseq > 10000) {
      print(cmd <- paste("vsearch -fastx_subsample", thisfile, "-sample_size 10000 -randseed 1234 -fastaout", outfile))
      system(cmd)
      nout <- 10000
    } else {
      file.copy(thisfile, outfile, overwrite = TRUE)
    }
    subsamptxt <- paste("Subsampled", nout, "of", nseq, "sequences")
    print(subsamptxt)
    # Save messages about number of subsampled sequences
    txtfile <- file.path(FASTAdir, paste0(thisRUNID, ".txt"))
    writeLines(subsamptxt, txtfile)
  }

  # Clean up
  file.remove(Sys.glob("*.fa"))
  file.remove(Sys.glob("*.fasta"))
  return()
}

# Find chimeras in combined FASTA file, then split result into samples 20200913
# (Saves time by masking and creating k-mer index only once)
findchimeras <- function() {
  message("=========================")
  message(paste0("Finding chimeras [", study, "]"))
  message("=========================")

  # Change to output directory
  olddir <- setwd(FASTAdir)
  on.exit(setwd(olddir))
  # Make combined FASTA file
  print(cmd <- "cat *.fa > combined.fasta")
  system(cmd)

  # Identify chimeras
  print(cmd <- paste("vsearch -uchime_ref combined.fasta -nonchimeras nonchimeras.fasta -db", refdb))
  system(cmd)

  # Split result into one file for each sample
  # https://stackoverflow.com/questions/26144692/printing-a-sequence-from-a-fasta-file
  allRUNID <- gsub("\\.fa$", "", dir(pattern = "\\.fa$"))
  for(thisRUNID in allRUNID) {
    # Don't use "awk '/>" because sample names for SBE+17 are at end of header of original FASTQ files, not beginning 20210501
    print(cmd <- paste0("awk '/", thisRUNID, "/{p++;print;next} /^>/{p=0} p' nonchimeras.fasta > ", thisRUNID, ".fasta"))
    system(cmd)
  }

  # Clean up
  file.remove("combined.fasta")
  file.remove("nonchimeras.fasta")
}

# Run RDP Classifier on one or multiple files 20200910
classify <- function(RUNID, conf = 0.8) {
  # Make sure RDP output directory exists 20200915
  if(!dir.exists(RDPdir)) stop(paste("directory", RDPdir, "does not exist"))
  # Force evaluation of RUNID before changing the directory
  allRUNID <- RUNID
  # Change to working directory
  olddir <- setwd(workdir)
  on.exit(setwd(olddir))

  # Loop over all RUNIDs
  n <- length(allRUNID)
  for(i in 1:n) {
    thisRUNID <- allRUNID[i]
    message("=======================================")
    message(paste0("Classifying ", thisRUNID, " [", study, "] (", i, "/", n, ")"))
    message("=======================================")

    # Input FASTA file
    FASTAfile <- file.path(FASTAdir, paste0(thisRUNID, ".fasta"))

    # Copy FASTA file to here without .fasta suffix (so it doesn't get into sample name)
    file.copy(FASTAfile, thisRUNID, overwrite = TRUE)
    # Run classifier with -c (confidence score threshold) 20201007
    print(cmd <- paste0("java -jar ", RDPjar, " classify -c ", conf, " -o ", thisRUNID, ".tab -h ", thisRUNID, ".txt ", thisRUNID))
    system(cmd)

    # Copy results (assignment count for each taxon) to output directories
    results <- paste0(thisRUNID, ".txt")
    resfile <- file.path(RDPdir, paste0(thisRUNID, ".txt"))
    file.copy(results, resfile, overwrite = TRUE)

    # Clean up
    file.remove(Sys.glob(paste0("*", thisRUNID, "*")))
  }
}

# Combine output files of RDP Classifier into a single CSV file for one study 20200903
# Use RDP Classifier merge-count command 20200910
# Get file list from RDPdir 20210817
mkRDP <- function() {
  # Change to temporary directory
  olddir <- setwd(tempdir())
  on.exit(setwd(olddir))
  # Remove any existing .txt files
  file.remove(Sys.glob("*.txt"))

  # Loop over runs
  runs <- gsub(".txt", "", dir(RDPdir))
  for(run in runs) {
    # The RDP result file
    RDPfile <- file.path(RDPdir, paste0(run, ".txt"))
    print(RDPfile)
    # Copy the RDP result file to here (temporary directory)
    file.copy(RDPfile, ".")
  }

  # Now list all the RDP result files
  files <- dir(pattern = "txt")
  # Name of output file
  outfile <- file.path(olddir, "RDP", paste0(study, ".tab"))
  # Remove output file in case it exists
  if(file.exists(outfile)) file.remove(outfile)
  # RDP jar file
  RDPjar <- "/opt/rdp_classifier_2.13/dist/classifier.jar"
  # Run merge-count command
  print(cmd <- paste("java -jar", RDPjar, "merge-count", outfile, paste(files, collapse = " ")))
  system(cmd)
}
