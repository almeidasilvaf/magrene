
# Load data and create example data sets
data(gma_paralogs)
data(gma_grn)
data(gma_ppi)

edgelist <- gma_grn[500:1000, 1:2]
paralogs_wgd <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]

# Start tests
test_that("find_lambda() returns motifs as a character vector", {
    motifs <- find_lambda(edgelist, paralogs_wgd)
    expect_equal(class(motifs), "character")
})

test_that("find_delta() returns motifs as a character vector", {
    motifs <- find_delta(edgelist, paralogs_wgd, gma_ppi)
    expect_true(is.null(motifs))
})

test_that("find_v() returns motifs as a character vector", {
    edgelist <- gma_grn[2000:4000, 1:2] # reducing for test purposes
    motifs <- find_v(edgelist, paralogs_wgd)
    expect_equal(class(motifs), "character")
})

test_that("find_bifan() returns motifs as an edge list", {
    paralogs <- rbind(
        paralogs_wgd,
        data.frame(duplicate1 = "Glyma.01G177200", 
                   duplicate2 = "Glyma.08G116700")
    )
    
    motifs <- find_bifan(gma_grn[600:1400, ], paralogs)
    expect_equal(class(motifs), "character")
})

