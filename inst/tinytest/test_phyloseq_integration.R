# Test added on 20230705

info <- "Different metrics can be chosen in ps_metrics()"
data(mouse.GTDB_207)
metrics <- c("GRAVY", "pI")
values <- ps_metrics(mouse.GTDB_207, metrics = metrics, refdb = "GTDB_207", quiet = TRUE)
expect_equal(colnames(values), metrics, info = info)

info <- "A single metric is still returned in a data frame"
metrics <- "Zc"
values <- ps_metrics(mouse.GTDB_207, metrics = metrics, refdb = "GTDB_207", quiet = TRUE)
expect_equal(colnames(values), metrics, info = info)
