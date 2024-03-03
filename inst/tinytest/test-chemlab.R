# Test added on 20240302
info <- "chemlab() returns an unknown variable, unchanged"
expect_equal(chemlab("xxx"), "xxx", info = info)
