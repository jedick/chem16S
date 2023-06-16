\docType{data}
\name{chem16S-data}
\alias{ps.silva}
\alias{ps.RDP}
\alias{ps.GTDB}
\title{phyloseq objects generated using the DADA2 Pipeline Tutorial}
\description{
  Objects in \code{\link{phyloseq-class}} format created by following the \href{https://benjjneb.github.io/dada2/tutorial.html}{DADA2 Pipeline Tutorial} (version 1.16, accessed on 2023-06-14).
}
\details{
  These files were made by using three different training sets to assign taxonomy to genus level:

  \describe{
    \item{\file{ps.silva}}{silva_nr99_v138.1: \href{https://zenodo.org/record/4587955}{Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2}}
    \item{\file{ps.RDP}}{rdp_train_set_18: \href{https://zenodo.org/record/4310151}{RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5)}}
    \item{\file{ps.GTDB}}{GTDB_bac120_arc53_ssu_r207: \href{https://zenodo.org/record/6655692}{GTDB: DADA2 formatted 16S rRNA gene sequences for both bacteria & archaea}}
  }
}
\examples{
# Load one data set
data(ps.silva)
# Summarize the object (phyloseq format)
ps.silva
}
\references{
Benjamin Callahan. DADA2 Pipeline Tutorial (1.16). \url{https://benjjneb.github.io/dada2/tutorial.html}, accessed on 2023-06-14.
}
\author{
Jeffrey M. Dick, following the DADA2 Pipeline Tutorial
}
\keyword{data}
