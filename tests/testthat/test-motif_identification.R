
# Load data and create example data sets
data(gma_paralogs)
data(gma_grn)
paralogs <- gma_paralogs[, 1:2]
gene1 <- "Glyma.01G011500"
gene2 <- "Glyma.01G011600"

edgelist <- gma_grn[1:500, 1:2]
paralogs_wgd <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]

# Start tests
test_that("are_paralogs() returns a logical scalar", {
    count <- count_lambda(edgelist, paralogs)
    expect_true(is.numeric(count))
})

