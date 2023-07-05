# chem16S/calc_metrics.R
# Calculate chemical metrics from amino acid compositions of proteins
# 20191027 initial version for canprot/metrics.R
# 20230704 adapted for chem16S

calc_metrics <- function(AAcomp, metrics = c("Zc", "nO2", "nH2O")) {

  values <- lapply(metrics, function(metric) {
  
    if(metric == "Zc") {

      # Calculate carbon oxidation state for amino acid compositions 20180228
      # The number of carbons of the amino acids
      nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
        Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
        Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
      # The Ztot of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula) * nC_AA
      Ztot_AA <- c(Ala = 0, Cys = 2, Asp = 4, Glu = 2, Phe = -4, Gly = 2, His = 4, 
        Ile = -6, Lys = -4, Leu = -6, Met = -2, Asn = 4, Pro = -2, Gln = 2, 
        Arg = 2, Ser = 2, Thr = 0, Val = -4, Trp = -2, Tyr = -2)
      # The Zc of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula)
      Zc_AA <- Ztot_AA / nC_AA
      # Find columns with names for the amino acids
      isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Zc_AA))
      iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(Zc_AA)))
      # Calculate the nC for all occurrences of each amino acid
      multC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
      # Multiply nC by Zc
      multZc <- t(t(multC) * Zc_AA[iAA])
      # Calculate the total Zc and nC, then the overall Zc
      Zctot <- rowSums(multZc)
      nCtot <- rowSums(multC)
      Zctot / nCtot

    } else if(metric == "nH2O") {

      # Calculate stoichiometric hydration state for proteins with given amino acid compositions 20181228
      basis <- "QEC"
      if(basis == "QEC") {
        # How to get the number of H2O in reactions to form amino acid residues from the "QEC" basis:
        ## library(CHNOSZ)
        ## basis("QEC")
        ## nH2O_AA <- species(aminoacids(""))$H2O
        ## names(nH2O_AA) <- aminoacids(3)
        nH2O_AA <- c( Ala =  0.6, Cys =    0, Asp = -0.2, Glu =    0, Phe = -2.2, Gly =  0.4, His = -1.8,
          Ile =  1.2, Lys =  1.2, Leu =  1.2, Met =  0.4, Asn = -0.2, Pro =    0, Gln =    0,
          Arg =  0.2, Ser =  0.6, Thr =  0.8, Val =    1, Trp = -3.8, Tyr = -2.2) - 1
        # Note: subtract 1 to get amino acid residues in proteins
      }
      # Find columns with names for the amino acids
      isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nH2O_AA))
      iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nH2O_AA)))
      # Calculate total number of H2O in reactions to form proteins
      nH2O <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * nH2O_AA[iAA]))
      ## Add one H2O for each pair of terminal groups (i.e., number of polypeptide chains)
      #terminal_H2O <- 0
      #nH2O <- nH2O + terminal_H2O
      # Divide by number of residues (length of protein)
      nH2O / rowSums(AAcomp[, isAA, drop = FALSE])

    } else if(metric == "nO2") {

      # Calculate stoichiometric oxidation state for proteins with given amino acid compositions 20201016
      basis <- "QEC"
      if(basis == "QEC") {
        # How to get the number of O2 in reactions to form amino acid residues from the "QEC" basis:
        ## library(CHNOSZ)
        ## basis("QEC")
        ## nO2_AA <- species(aminoacids(""))[["O2"]]
        ## names(nO2_AA) <- aminoacids(3)
        nO2_AA <- c(Ala = -0.3, Cys = 0, Asp = 0.6, Glu = 0, Phe = -1.9, Gly = 0.3, 
          His = 0.4, Ile = -2.1, Lys = -1.6, Leu = -2.1, Met = -1.2, Asn = 0.6, 
          Pro = -1, Gln = 0, Arg = -0.1, Ser = 0.2, Thr = -0.4, Val = -1.5, 
          Trp = -1.6, Tyr = -1.4)
      }
      # Find columns with names for the amino acids
      isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nO2_AA))
      iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nO2_AA)))
      # Calculate total number of O2 in reactions to form proteins
      nO2 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * nO2_AA[iAA]))
      # Divide by number of residues (length of protein)
      nO2 <- nO2 / rowSums(AAcomp[, isAA, drop = FALSE])

    } else if(metric == "GRAVY") {

      # Calculate GRAVY for amino acid compositions 20191024
      # Values of the hydropathy index from Kyte and Doolittle, 1982
      # doi:10.1016/0022-2836(82)90515-0
      Hind <- c(Ala =  1.8, Cys =  2.5, Asp = -3.5, Glu = -3.5, Phe =  2.8,
                Gly = -0.4, His = -3.2, Ile =  4.5, Lys = -3.9, Leu =  3.8,
                Met =  1.9, Asn = -3.5, Pro = -1.6, Gln = -3.5, Arg = -4.5,
                Ser = -0.8, Thr = -0.7, Val =  4.2, Trp = -0.9, Tyr = -1.3)
      # Find columns with names for the amino acids
      isAA <- colnames(AAcomp) %in% names(Hind)
      iAA <- match(colnames(AAcomp)[isAA], names(Hind))
      # Calculate total of hydropathy values for each protein
      sumHind <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Hind[iAA]))
      # Divide by length of proteins to get grand average of hydropathy (GRAVY)
      sumHind / rowSums(AAcomp[, isAA, drop = FALSE])

    } else if(metric == "pI") {

      # Calculate isoelectric point for proteins 20191026
      # A function to calculate isoelectric point for a single amino acid composition
      onepI <- function(AA) {
        # Find the column names of AAcomp that are in Ztab
        isZ <- names(AA) %in% dimnames(Ztab)[[2]]
        iZ <- match(names(AA)[isZ], dimnames(Ztab)[[2]])
        # Calculate the total charge as a function of pH
        # ... the "else" is in case we have a data frame (used when first writing this function)
        if(is.numeric(AA)) Ztot <- Ztab[, iZ] %*% AA[isZ]
        else Ztot <- Ztab[, iZ] %*% as.matrix(t(AA[, isZ]))
        # Find pH where charge is closest to zero
        # (absolute charge is minimized)
        ipH <- which.min(abs(Ztot))
        Ztab[ipH, 1]
      }
      # Number of N- and C-terminal groups is 1, unless the input data frame has a value for number of chains
      Nterm <- Cterm <- 1
      if(!is.null(AAcomp$chains)) Nterm <- Cterm <- AAcomp$chains
      if(!"Nterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Nterm = Nterm)
      if(!"Cterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Cterm = Cterm)
      # NOTE: apply() converts the input to matrix,
      # so we extract the numeric columns of AAcomp to avoid possible coercion of all values to character
      isnum <- unlist(lapply(AAcomp, "class")) %in% c("integer", "numeric")
      myAA <- AAcomp[, isnum, drop = FALSE]
      # Run the calculation for each composition
      apply(myAA, 1, onepI)

    } else if(metric == "MW") {

      # Calculate average molecular weight per amino acid residue 20200501
      # MW_AA <- sapply(CHNOSZ::makeup(info(aminoacids(""))), mass) - mass("H2O")
      # names(MW_AA) <- aminoacids(3)
      MW_AA <- c(Ala = 71.0788, Cys = 103.1388, Asp = 115.0886, Glu = 129.11548, 
        Phe = 147.17656, Gly = 57.05192, His = 137.14108, Ile = 113.15944, 
        Lys = 128.17408, Leu = 113.15944, Met = 131.19256, Asn = 114.10384, 
        Pro = 97.11668, Gln = 128.13072, Arg = 156.18748, Ser = 87.0782, 
        Thr = 101.10508, Val = 99.13256, Trp = 186.2132, Tyr = 163.17596
      )
      # Find columns with names for the amino acids
      isAA <- colnames(AAcomp) %in% names(MW_AA)
      iAA <- match(colnames(AAcomp)[isAA], names(MW_AA))
      # Calculate total MW of residues in each protein
      MW <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * MW_AA[iAA]))
      # Divide by number of residues (length of protein)
      MW / rowSums(AAcomp[, isAA, drop = FALSE])

    } else if(tolower(metric) == "length") {

      # Find columns with names for the amino acids
      AA_names <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
       "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
      isAA <- colnames(AAcomp) %in% AA_names
      # Protein length is number of amino acid residues 20230704
      length <- rowSums(AAcomp[, isAA, drop = FALSE])
      # Divide by number of proteins (polypeptide chains) if available 20230705
      if(!is.null(AAcomp$chains)) length <- length / AAcomp$chains
      length

    } else stop(paste0("'", metric, "' is not an available metric"))

  })

  values <- do.call(cbind, values)
  colnames(values) <- metrics
  as.data.frame(values)

}

#########################
### UNEXPORTED OBJECT ###
### ( used for pI )   ###
#########################

# Tabulate charges for sidechains and terminal groups from pH 0 to 14
Ztab <- local({
  # A function to calculate charge as a function of pH for a single group
  ZpH <- function(pK, Z, pH) {
    alpha <- 1/(1 + 10^(Z * (pH - pK)))
    alpha * Z
  }
  # List the pKs of the groups
  pK <- list(Cterm = 3.55, Nterm = 7.5,
    Asp = 4.05, Glu = 4.45, His = 5.98,
    Cys = 9, Lys = 10, Tyr = 10, Arg = 12
  )
  # List the unit charges of the groups
  Z <- list(Cterm = -1, Nterm = 1,
    Asp = -1, Glu = -1, His = 1,
    Cys = -1, Lys = 1, Tyr = -1, Arg = 1
  )
  # Get the charges for a range of pH values
  pH <- seq(0, 14, 0.01)
  Ztab <- mapply(ZpH, pK = pK, Z = Z, MoreArgs = list(pH = pH))
  # Add a column with the pH values
  cbind(pH = pH, Ztab)
})
