\encoding{UTF-8}
\name{calc_metrics}
\alias{calc_metrics}
\title{Calculate chemical metrics from amino acid composition}
\description{
Calculates various chemical metrics from amino acid composition(s) of one or more proteins.
}

\usage{
  calc_metrics(AAcomp, metrics = c("Zc", "nO2", "nH2O"))
}

\arguments{
  \item{AAcomp}{data frame of amino acid composition with column names being the 3-letter abbreviations of the amino acids (first letter uppercase)}
  \item{metrics}{character, chemical metrics to calculate}
}

\details{

The following metrics are available:
\itemize{
  \item{\code{Zc}:}{
    Carbon oxidation state.
    Its value is calculated from elemental composition as described by Dick (2014).
  }
  \item{\code{nO2}:}{
    Stoichiometric oxidation state.
    \code{nO2} and \code{nH2O} are calculated from theoretical formation reactions of proteins from the QEC basis species (glutamine, glutamic acid, cysteine, \H2O, and \O2) and are normalized by number of amino acid residues (see Dick et al., 2020).
  }
  \item{\code{nH2O}:}{
    Stoichiometric hydration state.
    See above; note that the contribution of terminal groups is not counted in the calculation of \code{nH2O}.
  }
  \item{\code{GRAVY}:}{
    Grand average of hydropathicity.
    Values of the hydropathy index for individual amino acids are from Kyte and Doolittle (1982).
  }
  \item{\code{pI}:}{
    Isoelectric point.
    The net charge for each ionizable group was pre-calculated from pH 0 to 14 at intervals of 0.01.
    The isoelectric point is found as the pH where the sum of charges of all groups in the protein is closest to zero.
    The \pK values for the terminal groups and sidechains are taken from Bjellqvist et al. (1993) and Bjellqvist et al. (1994); note that the calculation does not implement position-specific adjustments described in the latter paper.
    The number of N- and C-terminal groups is taken to be one, unless a column for \code{chains} (number of polypeptide chains) is given in \code{AAcomp}.
  }
  \item{\code{MW}:}{
    Molecular weight.
    This is the per-residue molecular weight and doesn't include the contribution of terminal groups.
  }
  \item{\code{length} or \code{Length}:}{
    Protein length.
    This is the number of amino acid residues per protein.
  }
  \item{\code{H/C}, \code{H_C}, or \code{HC}:}{
    H/C ratio (not counting terminal H-OH groups of the protein backbone).
  }
  \item{\code{N/C}, \code{N_C}, or \code{NC}:}{
    N/C ratio.
  }
  \item{\code{O/C}, \code{O_C}, or \code{OC}:}{
    O/C ratio (not counting terminal H-OH groups of the protein backbone).
  }
  \item{\code{S/C}, \code{S_C}, or \code{SC}:}{
    S/C ratio.
  }
}

}

\seealso{
\code{\link{get_metrics}}
}

\value{
A data frame with the same number of rows as \code{AAcomp} and one column of numeric values for each of the \code{metrics}.
An error is produced if any of the \code{metrics} is not available for calculation.
}

\examples{
# Amino acid composition of alanylglycine
AG <- data.frame(Ala = 1, Gly = 1)
# Calculate default metrics
calc_metrics(AG)
# Calculate selected metrics
calc_metrics(AG, c("GRAVY", "pI"))
}

\references{
Bjellqvist B, Hughes GJ, Pasquali C, Paquet N, Ravier F, Sanchez J-C, Frutiger S, Hochstrasser D. 1993. The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. \emph{Electrophoresis} \bold{14}: 1023--1031. \doi{10.1002/elps.11501401163}

Bjellqvist B, Basse B, Olsen E, Celis JE. 1994. Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a {pH} scale where isoelectric points correlate with polypeptide compositions. \emph{Electrophoresis} \bold{15}: 529--539. \doi{10.1002/elps.1150150171}

Dick JM. 2014. Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}: 20131095. \doi{10.1098/rsif.2013.1095}

Dick JM, Yu M, Tan J. 2020. Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}: 6145--6162. \doi{10.5194/bg-17-6145-2020}

Gasteiger E, Hoogland C, Gattiker A, Duvaud S, Wilkins MR, Appel RD, Bairoch A. 2005. Protein identification and analysis tools on the ExPASy server. In: Walker JM, editor. \emph{The Proteomics Protocols Handbook}. Totowa, NJ: Humana Press. pp. 571–607. \doi{10.1385/1-59259-890-0:571}

Kyte J, Doolittle RF. 1982. A simple method for displaying the hydropathic character of a protein. \emph{J. Mol. Biol.} \bold{157}: 105--132. \doi{10.1016/0022-2836(82)90515-0}
}
