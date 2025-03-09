# Test added on 20250309

# Clear the manual_mappings option
options(manual_mappings = NULL)
# Reset the option
chem16S:::.onLoad()
info <- ".onLoad() adds manual_mappings option"
mm <- getOption("manual_mappings")
expect_equal(dim(mm), c(16, 5))
