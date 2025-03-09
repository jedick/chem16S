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
expect_equal(p1$labels$x, "SampleType", info = info)
expect_equal(p1$labels$y, "Chemical Metric", info = info)
expect_equal(dim(p1$data), c(78, 10), info = info)
expect_equal(levels(p1$data$variable), c("Zc", "nO2", "nH2O"), info = info)

info <- "plot_ps_metrics2() returns plotting object with expected content"
p2 <- plot_ps_metrics2(GlobalPatterns, color = "SampleType", refdb = "RefSeq_206", quiet = TRUE)
expect_equal(names(p2$data)[1:2], c("Zc", "nH2O"), info = info)
expect_equal(dim(p2$data), c(26, 10), info = info)
