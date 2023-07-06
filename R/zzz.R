# chem16S/zzz.R

# Set 'manual_mappings' option 20230706
# Adapted from R/src/library/grDevices/zzz.R
.onLoad <- function(libname, pkgname) {
  manual_mappings <- read.csv(system.file("extdata/manual_mappings.csv", package = "chem16S"))
  op.chem16S <- list(manual_mappings = manual_mappings)
  toset <- !(names(op.chem16S) %in% names(.Options))
  if(any(toset)) options(op.chem16S[toset])
}
