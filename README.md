# chem16S
### Chemical metrics of estimated community proteomes from 16S rRNA data

This R package implements the method described in [Dick and Tan (2022)](https://doi.org/10.1007/s00248-022-01988-9) for combining **RDP Classifier** output with **reference proteomes** of archaea and bacteria to generate the amino acid compositions of **estimated community proteomes**.

The amino acid compositions of the estimated community proteomes can be used to calculate chemical metrics such as **carbon oxidation state** (*Z*<sub>C</sub>) and **stoichiometric hydration state** (*n*<sub>H<sub>2</sub>O</sub>).
This example from the [help page for the `plotmet` function](man/plotmet.Rd) shows relatively low *n*<sub>H<sub>2</sub>O</sub> in high-salinity samples from the Baltic Sea obtained by processing 16S rRNA gene sequencing data from [Herlemann et al. (2016)](https://doi.org/10.3389/fmicb.2016.01883).

<!-- Default image is too big
![plotmet example: Baltic Sea nH2O-ZC plot](inst/images/plotmet.png)
-->
<img src="inst/images/plotmet.png" width="500" />

### Installation

Install the **remotes** package from CRAN, then use that to install **chem16S** from GitHub.

```
install.packages("remotes")
remotes::install_github("jedick/chem16S")
```
