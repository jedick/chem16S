#/bin/sh
# gencat.sh: extract microbial protein records from the RefSeq catalog
# 20130120 First(?) version for RefSeq release 57 (in CHNOSZ)
# 20190830 Update for RefSeq release 95 (use accessions not gi numbers)
# 20200716 Update for RefSeq release 201
# 20210523 Version number update for RefSeq release 206

# zcat: Decompress the release catalog
# awk: Save the tax_id, accession, and FTP directory columns
#      (We must drop the taxon name because e.g. Staphylococcus aureus FP_N5208 gives a false match in grep)
# grep: Extract the microbial protein records (169527530 records, 6.3 GB)
zcat RefSeq-release206.catalog.gz | awk -F$'\t' '{print $1,$3,$4}' | grep "P_.*bacteria\|P_.*archaea\|P_.*viral" > microbial.protein.tmp
# Save only the accession and taxid columns (3.5 GB)
awk '{print $2,$1}' microbial.protein.tmp > accession.taxid.tmp
rm microbial.protein.tmp
# Sort the file on accession so that it can be used with the unix 'join' command
# NOTE: for join to work, -n for numeric sort is *not* used here
sort accession.taxid.tmp > accession.taxid.txt
rm accession.taxid.tmp

# Make a list of unique taxids 20190830
# For RefSeq 206, the number is 49298
awk '{print $2}' accession.taxid.txt > taxid.tmp
sort -u taxid.tmp > unique_taxids.txt
# Also list the unique WP accessions
# For RefSeq 206, the number is 33584
awk '{print $2}' release206.MultispeciesAutonomousProtein2taxname > taxid.tmp
sort -u taxid.tmp > WP_unique_taxids.txt
rm taxid.tmp
