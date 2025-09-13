# Tests added on 20230705

info <- "Different metrics can be chosen in ps_metrics()"
data(mouse.GTDB_220)
metrics <- c("GRAVY", "pI")
values <- ps_metrics(mouse.GTDB_220, metrics = metrics, refdb = "GTDB_220", quiet = TRUE)
expect_equal(colnames(values), metrics, info = info)

info <- "A single metric is still returned in a data frame"
metrics <- "Zc"
values <- ps_metrics(mouse.GTDB_220, metrics = metrics, refdb = "GTDB_220", quiet = TRUE)
expect_equal(colnames(values), metrics, info = info)

# Tests added on 20250309

# Start NULL graphics device
pdf(NULL)

info <- "plot_ps_metrics() returns plotting object with expected content"
data(GlobalPatterns, package = "phyloseq")
p1 <- plot_ps_metrics(GlobalPatterns, x = "SampleType", sortby = "Zc", refdb = "RefSeq_206")
# Build the plot object to get labels 20250913
# Needed after https://github.com/tidyverse/ggplot2/pull/5879
p1_b <- ggplot2::ggplot_build(p1)
expect_equal(p1_b$plot$labels$x, "SampleType", info = info)
expect_equal(p1_b$plot$labels$y, "Chemical Metric", info = info)
expect_equal(dim(p1$data), c(78, 10), info = info)
expect_equal(levels(p1$data$variable), c("Zc", "nO2", "nH2O"), info = info)

info <- "plot_ps_metrics2() returns plotting object with expected content"
p2 <- plot_ps_metrics2(GlobalPatterns, color = "SampleType", refdb = "RefSeq_206", quiet = TRUE)
expect_equal(names(p2$data)[1:2], c("Zc", "nH2O"), info = info)
expect_equal(dim(p2$data), c(26, 10), info = info)

# 20250321
info <- "plot_ps_metrics() handles 'color' and 'shape' arguments"
p3 <- plot_ps_metrics(GlobalPatterns, x = "SampleType", color = "SampleType", sortby = "Zc", refdb = "RefSeq_206")
p3_b <- ggplot2::ggplot_build(p3)
expect_equal(p3_b$plot$labels$colour, "SampleType", info = info)
p4 <- plot_ps_metrics(GlobalPatterns, x = "SampleType", shape = "SampleType", sortby = "Zc", refdb = "RefSeq_206")
p4_b <- ggplot2::ggplot_build(p4)
expect_equal(p4_b$plot$labels$shape, "SampleType", info = info)

info <- "plot_ps_metrics2() stops for too few metrics and warns for too many"
expect_error(plot_ps_metrics2(GlobalPatterns, metrics = c("Zc"), color = "SampleType", refdb = "RefSeq_206", quiet = TRUE),
  "please supply the names of two metrics", info = info)
expect_warning(plot_ps_metrics2(GlobalPatterns, metrics = c("Zc", "nH2O", "nO2"), color = "SampleType", refdb = "RefSeq_206", quiet = TRUE),
  "using the first two", info = info)

dev.off()
