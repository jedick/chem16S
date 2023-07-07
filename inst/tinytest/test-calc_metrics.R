## Tests for Zc, nO2, and nH2O adapted from canprot/tests/test-metrics.R on 20230704

info <- "Results are as expected for Zc, nO2, and nH2O"

## Calculate metrics for a few proteins the "long way" (using functions in CHNOSZ)
#library(CHNOSZ)
#basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
#Zc.ref <- ZC(protein.formula(1:6))
#nO2.ref <- protein.basis(1:6)[, "O2"] / protein.length(1:6)
## NOTE: subtract 1 so as exclude terminal groups from calculation of nH2O
#nH2O.ref <- (protein.basis(1:6)[, "H2O"] - 1) / protein.length(1:6)

Zc.ref <- c(-0.11633875106929, -0.0272787757817698, -0.195689166193988, -0.0492957746478873, -0.170212765957447, 0.0163132137030995)
nO2.ref <- c(-0.699539170506912, -0.522294022617124, -0.81466049382716, -0.574137931034483, -0.716346153846154, -0.471317829457364)
nH2O.ref <- c(-1.17465437788018, -0.881098546042003, -0.941666666666667, -0.955172413793103, -0.730769230769231, -0.886821705426357)

# Calculate metrics using calc_metrics() function in chem16S
AAcomp <- 
structure(list(protein = c("O08452", "AMY", "AMYA", "BPT1", "CYC", 
"LYSC"), organism = c("PYRFU", "BACSU", "PYRFU", "BOVIN", "BOVIN", 
"CHICK"), ref = c("UniProt", "UniProt", "UniProt", "UniProt", 
"UniProt", "UniProt"), abbrv = c("O08452", "P00691", "P49067", 
"P00974", "P62894", "P00698"), chains = c(1L, 1L, 1L, 1L, 1L, 
1L), Ala = c(28, 49, 26, 6, 6, 12), Cys = c(5, 1, 2, 6, 2, 8), 
    Asp = c(33, 44, 35, 2, 3, 7), Glu = c(23, 23, 66, 2, 9, 2
    ), Phe = c(20, 20, 37, 4, 4, 3), Gly = c(45, 51, 44, 6, 14, 
    12), His = c(12, 16, 14, 0, 3, 1), Ile = c(25, 35, 41, 2, 
    6, 6), Lys = c(19, 30, 48, 4, 18, 6), Leu = c(27, 36, 59, 
    2, 6, 8), Met = c(4, 10, 12, 1, 2, 2), Asn = c(21, 54, 24, 
    3, 5, 14), Pro = c(20, 23, 28, 4, 4, 2), Gln = c(7, 29, 15, 
    1, 3, 3), Arg = c(14, 24, 35, 6, 2, 11), Ser = c(21, 55, 
    33, 1, 1, 10), Thr = c(16, 45, 12, 3, 8, 7), Val = c(31, 
    32, 59, 1, 3, 6), Trp = c(26, 14, 17, 0, 1, 6), Tyr = c(37, 
    28, 41, 4, 4, 3)), row.names = c(NA, 6L), class = "data.frame")

metrics <- calc_metrics(AAcomp)

# Perform the tests
expect_equivalent(metrics$Zc, Zc.ref, info = info)
expect_equivalent(metrics$nO2, nO2.ref, info = info)
expect_equivalent(metrics$nH2O, nH2O.ref, info = info)

## Tests for GRAVY, pI, and MW adapted from canprot/man/metrics.Rd on 20230704

info <- "Results are as expected for GRAVY"

# Reference values obtained with ProtParam on uniprot.org
# https://web.expasy.org/cgi-bin/protparam/protparam1?P00698@19-147@
# https://web.expasy.org/cgi-bin/protparam/protparam1?P61823@27-150@
# https://web.expasy.org/cgi-bin/protparam/protparam1?P49067@2-649@
GRAVY.ref <- c(-0.472, -0.663, -0.325)

## Get the protein index in CHNOSZ's list of proteins
#iprotein <- pinfo(c("LYSC_CHICK", "RNAS1_BOVIN", "AMYA_PYRFU", "CSG_HALJP"))
## Then get the amino acid compositions
#AAcomp <- pinfo(iprotein)

AAcomp <-
structure(list(protein = c("LYSC", "RNAS1", "AMYA", "CSG"), organism = c("CHICK", 
"BOVIN", "PYRFU", "HALJP"), ref = c("UniProt", "UniProt", "UniProt", 
"UniProt"), abbrv = c("P00698", "P61823", "P49067", "Q9C4B4"), 
    chains = c(1L, 1L, 1L, 1L), Ala = c(12, 12, 26, 61), Cys = c(8, 
    8, 2, 0), Asp = c(7, 5, 35, 122), Glu = c(2, 5, 66, 86), 
    Phe = c(3, 3, 37, 20), Gly = c(12, 3, 44, 78), His = c(1, 
    4, 14, 4), Ile = c(6, 3, 41, 47), Lys = c(6, 10, 48, 4), 
    Leu = c(8, 2, 59, 47), Met = c(2, 4, 12, 0), Asn = c(14, 
    10, 24, 51), Pro = c(2, 4, 28, 29), Gln = c(3, 7, 15, 21), 
    Arg = c(11, 4, 35, 19), Ser = c(10, 15, 33, 74), Thr = c(7, 
    10, 12, 79), Val = c(6, 9, 59, 66), Trp = c(6, 0, 17, 2), 
    Tyr = c(3, 6, 41, 18)), row.names = c(6L, 9L, 3L, 14L), class = "data.frame")

# Calculate GRAVY
metrics <- calc_metrics(AAcomp[1:3, ], "GRAVY")
# Perform the test
expect_equivalent(round(metrics$GRAVY, 3), GRAVY.ref, info = info)

info <- "Results are as expected for pI"

# Reference values calculated with ProtParam on uniprot.org
# LYSC_CHICK: residues 19-147 (sequence v1)
# RNAS1_BOVIN: residues 27-150 (sequence v1)
# AMYA_PYRFU: residues 2-649 (sequence v2)
# CSG_HALJP: residues 35-862 (sequence v1)
pI.ref <- c(9.32, 8.64, 5.46, 3.37)

# Calculate pI
metrics <- calc_metrics(AAcomp, "pI")
# Perform the test
expect_equivalent(metrics$pI, pI.ref, info = info)

info <- "Results are as expected for MW and length"

MW.ref <- c(14313.14, 13690.29, 76178.25)
metrics <- calc_metrics(AAcomp, c("MW", "length"))
# MW gives the molecular weight per residue;
# multiply by length to get molecular weight per protein
MW.calc <- metrics$MW * metrics$length
# Add molecular weight of terminal H- and -OH groups
MW.calc <- MW.calc + 18.01528
expect_equivalent(round(MW.calc[1:3], 2), MW.ref, info = info)

## Tests for H/C, N/C, O/C, and S/C added on 20230707

# library(CHNOSZ)
# pf <- as.data.frame(protein.formula(AAcomp))
# pf$H <- pf$H - 2  # Remove terminal H-OH 
# pf$O <- pf$O - 1  # Remove terminal H-OH 
# HCref <- pf$H / pf$C
# OCref <- pf$O / pf$C
# NCref <- pf$N / pf$C
# SCref <- pf$S / pf$C

HC.ref <- c(1.56117455138662, 1.57739130434783, 1.50964265456608, 1.53856636685745)
OC.ref <- c(0.300163132137031, 0.333913043478261, 0.276517300056721, 0.405287544289997)
NC.ref <- c(0.314845024469821, 0.297391304347826, 0.250992626205332, 0.264649768329245)
SC.ref <- c(0.0163132137030995, 0.0208695652173913, 0.00397050482132728, 0)

metrics <- calc_metrics(AAcomp, c("HC", "OC", "NC", "SC"))
expect_equivalent(metrics$HC, HC.ref)
expect_equivalent(metrics$NC, NC.ref)
expect_equivalent(metrics$OC, OC.ref)
expect_equivalent(metrics$SC, SC.ref)
