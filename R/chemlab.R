# chem16S/chemlab.R
# Function to format figure labels
# Adapted from canprot's cplab.R on 20230704
# Wrap cplab rather than copying the code
#   (DRY - Don't Repeat Yourself) 20240302

chemlab <- function(varname) {
  out <- cplab[[varname]]
  if(is.null(out)) out <- varname
  out
}
