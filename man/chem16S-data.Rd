\docType{data}
\name{chem16S-data}
\alias{mouse.silva}
\alias{mouse.RDP}
\alias{mouse.GTDB_214}
\title{`phyloseq-class` objects generated using the DADA2 Pipeline Tutorial}
\description{
  Objects in \code{\link[phyloseq]{phyloseq-class}} format created by following the DADA2 Pipeline Tutorial (Callahan, 2020).
}

\details{
  The example dataset for gut communities in a single mouse was taken from the \href{https://mothur.org/wiki/miseq_sop/}{mothur MiSeq SOP}.
  It is an extract of the complete dataset reported by Schloss et al. (2012).
  The output files were generated by using three different training sets to assign taxonomy to genus level:

  \describe{
    \item{\file{mouse.silva}}{silva_nr99_v138.1: \href{https://zenodo.org/record/4587955}{Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2}}
    \item{\file{mouse.RDP}}{rdp_train_set_18: \href{https://zenodo.org/record/4310151}{RDP taxonomic training data formatted for DADA2 (RDP trainset 18/release 11.5)}}
    \item{\file{mouse.GTDB_214}}{GTDB_bac120_arc53_ssu_r214: \href{https://zenodo.org/record/10403693}{GTDB: DADA2 formatted 16S rRNA gene sequences for both bacteria & archaea}}
  }
}

\examples{
# Load one data set
data(mouse.silva)
# Summarize the object (phyloseq-class)
mouse.silva
}

\references{
Callahan B (2020). DADA2 Pipeline Tutorial (1.16). \url{https://benjjneb.github.io/dada2/tutorial.html}, accessed on 2023-06-14. Publication year taken from \url{https://code.bioconductor.org/browse/dada2/blob/RELEASE_3_11/DESCRIPTION}.

Schloss PD, Schubert AM, Zackular JP, Iverson KD, Young VB, Petrosino JF. 2012. Stabilization of the murine gut microbiome following weaning. \emph{Gut Microbes} \bold{3}(4): 383–393. \doi{10.4161/gmic.21008}
}

\author{
Jeffrey M. Dick, following the DADA2 Pipeline Tutorial by Callahan (2020).
}

\keyword{data}
