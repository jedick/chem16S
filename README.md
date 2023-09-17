<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/chem16S)](https://cran.r-project.org/package=chem16S)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793059.svg)](https://doi.org/10.5281/zenodo.6793059)
[![R-CMD-check](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

#### *chem16S* processes taxonomic classifications of 16S rRNA gene sequences for a chemical representation of the genomic differences between communities.

* Read the paper in *Bioinformatics*: [*chem16S*: Community-level chemical metrics for exploring genomic adaptation to environments](https://doi.org/10.1093/bioinformatics/btad564).

* View the [manual](https://chnosz.net/chem16S/manual) and vignettes: [Chemical metrics of reference proteomes for taxa](https://chnosz.net/chem16S/vignettes/metrics.html), [Integration of *chem16S* with *phyloseq*](https://chnosz.net/chem16S/vignettes/phyloseq.html), and [Plotting two chemical metrics](https://chnosz.net/chem16S/vignettes/plotting.html).

Supported input formats:
* `phyloseq-class` objects created using [*phyloseq*](https://doi.org/doi:10.18129/B9.bioc.phyloseq)
* [RDP Classifier](https://sourceforge.net/projects/rdp-classifier/)

Supported reference databases:

* [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) (GTDB release 207)
* [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/) (RefSeq release 206)

The *chem16S* R package combines taxonomic classifications of high-throughput 16S rRNA gene sequences with precomputed amino acid compositions of reference proteomes for archaea and bacteria to generate the amino acid compositions of **community reference proteomes**.
Chemical metrics of community reference proteomes such as **carbon oxidation state** (*Z*<sub>C</sub>) and **stoichiometric hydration state** (*n*<sub>H<sub>2</sub>O</sub>) reveal new ways that microbial genomes are adapted to environmental conditions.
For instance, an association of lower *n*<sub>H<sub>2</sub>O</sub> with higher salinity in the Baltic Sea suggests a genomically encoded dehydration trend:

<!-- Default image is too big
![chem16S::plot_metrics example: Baltic Sea nH2O-Zc plot](inst/images/plot_metrics.png)
-->
<img src="inst/images/plot_metrics.png" alt="Baltic Sea nH2O-Zc plot (example from chem16S::plot_metrics)" width="400" />

PSU stands for practical salinity units.
The sequence data analyzed for this plot was taken from [Herlemann et al. (2016)](https://doi.org/10.3389/fmicb.2016.01883) and the code to make this plot is available in the [help page for `chem16S::plot_metrics`](https://chnosz.net/chem16S/manual/plot_metrics.html).

### Methods

* Scripts in the [GTDB](inst/extdata/GTDB) and [RefSeq](inst/extdata/RefSeq) directories were used to generate the reference proteomes for genus- and higher-level archaeal and bacterial taxa (and viruses for RefSeq).

* It is recommended to use a GTDB training set for taxonomic classification (e.g. [*DADA2 formatted 16S rRNA gene sequences for both bacteria & archaea*](https://doi.org/10.5281/zenodo.6655692)) so that taxonomic assignments can be automatically matched to GTDB reference proteomes available in *chem16S*.

* For taxonomic classifications made using the RDP training set, *chem16S* includes manual mappings to the NCBI taxonomy described by [Dick and Tan (2023)](https://doi.org/10.1007/s00248-022-01988-9).

### Installation

First install *phyloseq* from Bioconductor:

```
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
```

Then install the release version of *chem16S* from CRAN:

```
install.packages("chem16S")
```

Or use `install_github` from *remotes* or *devtools* to install the development version of *chem16S* from GitHub:

```
if(!require("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jedick/chem16S", build_vignettes = TRUE)
```
