<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793059.svg)](https://doi.org/10.5281/zenodo.6793059)
[![R-CMD-check](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# chem16S

### Chemical metrics of 16S rRNA-based community reference proteomes

This R package combines **taxonomic classifications** of high-throughput 16S rRNA gene sequences with **reference proteomes** of archaea and bacteria to generate the amino acid compositions of **community reference proteomes**. Taxonomic classifications can be read from the output of the [RDP Classifier](https://sourceforge.net/projects/rdp-classifier/) or from `phyloseq-class` objects created using the Bioconductor package [**phyloseq**](https://doi.org/doi:10.18129/B9.bioc.phyloseq).

The amino acid compositions of community reference proteomes are used to calculate chemical metrics such as **carbon oxidation state** (*Z*<sub>C</sub>) and **stoichiometric hydration state** (*n*<sub>H<sub>2</sub>O</sub>).
Lower *n*<sub>H<sub>2</sub>O</sub> is associated with increasing salinity in samples from the Baltic Sea:

<!-- Default image is too big
![chem16S::plot_metrics example: Baltic Sea nH2O-Zc plot](inst/images/plot_metrics.png)
-->
<img src="inst/images/plot_metrics.png" alt="chem16S::plot_metrics example: Baltic Sea nH2O-Zc plot" width="500" />

The code to make this plot is from the [help page for the `plot_metrics` function](man/plot_metrics.Rd) and uses sequence data reported by [Herlemann et al. (2016)](https://doi.org/10.3389/fmicb.2016.01883).

### Reference proteomes for taxa

Precomputed amino acid compositions of reference proteomes are provided for the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) (GTDB release 207) and the [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/) (RefSeq release 206). See the files in [`inst/extdata/RefSeq`](inst/extdata/RefSeq) for the steps used to download protein sequences from RefSeq and calculate the total amino acid composition for each NCBI taxonomic ID (taxid). The `taxon_AA.R` scripts for [GTDB](inst/extdata/GTDB) and [RefSeq](inst/extdata/RefSeq) were used to generate the reference proteomes for genus- and higher-level archaeal and bacterial taxa (and viruses for RefSeq).

* For 16S rRNA gene sequences classified with a GTDB-based training set (see [*DADA2 formatted 16S rRNA gene sequences for both bacteria & archaea*](https://doi.org/10.5281/zenodo.6655692)), taxonomic ranks and names are matched to the GTDB taxonomy to look up reference proteomes.

* For 16S rRNA gene sequences classified using the default training set of the RDP Classifier, **chem16S** includes manual mapping to the RefSeq taxonomy as described by [Dick and Tan (2023)](https://doi.org/10.1007/s00248-022-01988-9).

### Installation

After installing **phyloseq** from Bioconductor, use `install_github` (provided by either **remotes** or **devtools**) to install **chem16S** from GitHub.

```
# Install 'BiocManager' from CRAN
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install 'phyloseq' from Bioconductor
BiocManager::install("phyloseq")

# Install 'remotes' from CRAN
if(!require("remotes", quietly = TRUE)) install.packages("remotes")

# Install 'chem16S' from GitHub
remotes::install_github("jedick/chem16S", build_vignettes = TRUE)
```
