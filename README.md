[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793059.svg)](https://doi.org/10.5281/zenodo.6793059)

# chem16S

### Chemical metrics of community reference proteomes from 16S rRNA data

This R package implements the method described in [Dick and Tan (2023)](https://doi.org/10.1007/s00248-022-01988-9) for combining **RDP Classifier** output with **reference proteomes** of archaea and bacteria to generate the amino acid compositions of **community reference proteomes**.

The amino acid compositions of the community reference proteomes can be used to calculate chemical metrics such as **carbon oxidation state** (*Z*<sub>C</sub>) and **stoichiometric hydration state** (*n*<sub>H<sub>2</sub>O</sub>).
This example from the help page for the [`plotmet`](man/plotmet.Rd) function shows relatively low *n*<sub>H<sub>2</sub>O</sub> in high-salinity samples from the Baltic Sea with 16S rRNA gene sequencing data from [Herlemann et al. (2016)](https://doi.org/10.3389/fmicb.2016.01883).

<!-- Default image is too big
![chem16S::plotmet example: Baltic Sea nH2O-ZC plot](inst/images/plotmet.png)
-->
<img src="inst/images/plotmet.png" alt="chem16S::plotmet example: Baltic Sea nH2O-ZC plot" width="500" />

### Reference proteomes for taxa

See README.txt and scripts in [`inst/extdata/refseq`](inst/extdata/refseq) for steps used to download protein sequences from the RefSeq database and calculate the total amino acid composition for each NCBI taxonomic ID (taxid).

The `taxon_AA.R` script in [`inst/extdata/chem16S`](inst/extdata/chem16S) generates the reference proteomes for each viral, archaeal and bacterial genus, family, order, class, and phylum in the RefSeq database as follows:

* Only taxids classified at the species level are used, and archaeal and bacterial species with less than 500 reference protein sequences are excluded;
* For each species-level taxid, the total amino acid composition is converted to per-protein mean amino acid composition (this is done so that species with different proteome sizes contribute equally to the reference proteomes of higher-level taxa);
* For each genus, the mean amino acid compositions of all species-level taxids in that genus are summed and divided by the number of taxids to get the amino acid composition of the reference proteome;
* Analogously, the mean amino acid compositions of all species-level taxids in each family, order, class, and phylum are used to get the reference proteomes for taxa at those levels.

### Installation

Install the **remotes** package from CRAN, then use that to install **chem16S** from GitHub.

```
install.packages("remotes")
remotes::install_github("jedick/chem16S")
```
