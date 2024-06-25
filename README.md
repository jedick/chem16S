<!-- badges: start -->
[![CRAN](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fcloud.r-project.org%2Fweb%2Fpackages%2Fchem16S%2FDESCRIPTION&query=%24.Version&logo=r&label=CRAN&color=4bc51e)](https://cran.r-project.org/package=chem16S)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793059.svg)](https://doi.org/10.5281/zenodo.6793059)
[![R-CMD-check](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jedick/chem16S/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

*chem16S* calculates chemical metrics for microbial communities by combining taxonomic abundances with genomic reference sequences for proteins.
The chemical representation of communities has applications ranging from human microbiomes to Earth-life coevolution.

* Read the paper in *Bioinformatics*: [*chem16S*: community-level chemical metrics for exploring genomic adaptation to environments](https://doi.org/10.1093/bioinformatics/btad564).

* View the [manual](https://chnosz.net/chem16S/manual/) and vignettes: [Chemical metrics of reference proteomes for taxa](https://chnosz.net/chem16S/vignettes/metrics.html), [Integration of *chem16S* with *phyloseq*](https://chnosz.net/chem16S/vignettes/phyloseq.html), and [Plotting two chemical metrics](https://chnosz.net/chem16S/vignettes/plotting.html).

Supported input formats:
* `phyloseq-class` objects created using [*phyloseq*](https://doi.org/doi:10.18129/B9.bioc.phyloseq)
* [RDP Classifier](https://sourceforge.net/projects/rdp-classifier/)

Supported reference databases:

* [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) (GTDB release 214)
* [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/) (RefSeq release 206)

### Description

The *chem16S* R package combines taxonomic classifications of high-throughput 16S rRNA gene sequences with precomputed amino acid compositions of reference proteomes for archaea and bacteria to obtain the amino acid compositions of **community reference proteomes**.
Chemical metrics of community reference proteomes such as **carbon oxidation state** (*Z*<sub>C</sub>) and **stoichiometric hydration state** (*n*<sub>H<sub>2</sub>O</sub>) reveal new types of adaptations of microbial genomes to environmental conditions.
For instance, an association of lower *n*<sub>H<sub>2</sub>O</sub> with higher salinity in the Baltic Sea suggests a genomically encoded dehydration trend:

<!-- Default image is too big
![chem16S::plot_metrics example: Baltic Sea nH2O-Zc plot](inst/images/plot_metrics.png)
-->
<a href="https://chnosz.net/chem16S/manual/plot_metrics.html"><img src="inst/images/plot_metrics.png" alt="Baltic Sea nH2O-Zc plot (example from chem16S::plot_metrics)" width="400" /></a>

PSU stands for practical salinity units.
The sequence data analyzed for this plot was taken from [Herlemann et al. (2016)](https://doi.org/10.3389/fmicb.2016.01883) and the code to make this plot is available in the [help page for `chem16S::plot_metrics`](https://chnosz.net/chem16S/manual/plot_metrics.html).

### Methods

* Scripts in the [GTDB_214](inst/RefDB/GTDB_214) and [RefSeq_206](inst/RefDB/RefSeq_206) directories were used to generate reference proteomes for genus- and higher-level archaeal and bacterial taxa (and viruses for RefSeq).

* It is recommended to use *16S rRNA sequences from GTDB* for taxonomic classification (files are available for [DADA2](https://doi.org/10.5281/zenodo.10403693) and the [RDP Classifier](https://doi.org/10.5281/zenodo.12525163)) so that taxonomic assignments can be automatically matched to GTDB reference proteomes available in *chem16S*.

* For taxonomic classifications made using the RDP training set (No. 18 07/2020, used in RDP Classifier version 2.13), *chem16S* includes manual mappings to the NCBI taxonomy described by [Dick and Tan (2023)](https://doi.org/10.1007/s00248-022-01988-9).

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
