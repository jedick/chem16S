# Run R scripts to produce output files in this order:
# 1. genome_AA.R --> genome_AA.csv (hosted in https://github.com/jedick/JMDplots because of its size)
# 2. taxonomy.R  --> taxonomy.csv (mirrored in https://github.com/jedick/JMDplots)
# 3. taxon_AA.R  --> taxon_AA.csv

### Details for steps 1 and 2 ###

# This document describes steps for processing the 
# RefSeq database (release 206, 2021-05-21) to produce these files:

genome_AA.csv: total amino acid composition of all proteins for
  bacteria, archaea, and viral genomes in the RefSeq collection (n = 49448)
taxonomy.csv: taxid, phylum name and species name for 49448 taxa

# These functions/scripts have the following purpose (output files listed in parentheses):
gencat.sh - extract accession, taxid, sequence length from RefSeq release catalog (accession.taxid.txt)
genome_AA.R - get average amino acid composition for each taxid from gzipped sequence files (genome_AA.csv)
taxonomy.R - get taxonomic names for each taxid represented (taxonomy.csv)

## Download stuff

1. Download catalog files
wget https://ftp.ncbi.nih.gov/refseq/release/release-catalog/RefSeq-release206.catalog.gz (1.5 GB)
wget https://ftp.ncbi.nih.gov/refseq/release/release-catalog/release206.MultispeciesAutonomousProtein2taxname.gz (229 MB)
wget https://ftp.ncbi.nih.gov/refseq/release/release-catalog/release206.files.installed (1.5 MB)

2. List URLS for the microbial protein sequence files
grep bacteria.*faa release206.files.installed | awk '!($1="")' | sed -e "s/^\ /https\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/bacteria\//g" > urllist
grep archaea.*faa release206.files.installed | awk '!($1="")' | sed -e "s/^\ /https\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/archaea\//g" >> urllist
grep viral.*faa release206.files.installed | awk '!($1="")' | sed -e "s/^\ /https\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/viral\//g" >> urllist

3. Download the sequence files; use -N to skip existing files [1297 files, 33 GB]
wget -N -i urllist

4. Move the *.faa.gz files to a directory named 'protein'

## Protein stuff

# Timings were made on a 2018 Intel i7 laptop
5. gzip -d release206.MultispeciesAutonomousProtein2taxname.gz [2.8 GB, needed for finding all taxids for WP accessions]

6. Use 'gencat.sh' to generate accession.taxid.txt for bacteria, archaea, and viral proteins in the catalog [11 minutes]
   for RefSeq 206, 'wc -l accession.taxid.txt' is 169527530

7. Generate genome_AA.csv
   mkdir csv
   Then, run these R commands in the working directory that contains accession.taxid.txt, protein/, and csv/
   > source("genome_AA.R")
   > read.allfiles()    # ca. 19 hours (8 cores)
   > genome_AA()        # 5 minutes

## Taxonomy stuff

8. Edit 'taxonomy.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. These files can be downloaded from
    https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

    File info on server (accessed on 2021-05-21):
    taxdump.tar.gz            2021-05-21 21:27   53M

9. Source 'taxonomy.R' to generate the file 'taxonomy.csv' [8 hours]
