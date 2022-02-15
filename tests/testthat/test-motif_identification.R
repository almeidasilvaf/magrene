
# Load data and create example data sets
data(gma_paralogs)
data(gma_grn)
data(gma_ppi)

paralogs <- gma_paralogs[, 1:2]
gene1 <- "Glyma.01G011500"
gene2 <- "Glyma.01G011600"

edgelist <- gma_grn[500:1000, 1:2]
paralogs_wgd <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]

# Start tests
test_that("are_paralogs() returns a logical scalar", {
    p <- are_paralogs(gene1, gene2, paralogs)
    expect_true(is.logical(p))
})

test_that("are_interacting() returns a logical scalar", {
    i <- are_interacting("Glyma.19G213200", "Glyma.01G004300", gma_ppi)
    expect_true(is.logical(i))
})

test_that("find_lambda() returns motifs as an edge list", {
    motifs <- find_lambda(edgelist, paralogs)
    expect_equal(class(motifs), "list")
})

test_that("find_delta() returns motifs as an edge list", {
    motifs <- find_delta(edgelist, paralogs, gma_ppi)
    expect_true(is.null(motifs))
})
