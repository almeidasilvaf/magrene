
# Load data and create example data sets
data(gma_paralogs)
data(gma_grn)
data(gma_ppi)

edgelist <- gma_grn[500:1000, 1:2]
paralogs_wgd <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]

# Start tests
test_that("find_lambda() returns motifs as a character vector", {
    motifs <- find_lambda(edgelist, paralogs_wgd)
    motifs_count <- find_lambda(edgelist, paralogs_wgd, count_only = TRUE)
    
    expect_equal(class(motifs), "character")
    expect_true(is.numeric(motifs_count))
    expect_equal(length(motifs_count), 1)
})

test_that("find_delta() returns motifs as a character vector", {
    motifs <- find_delta(edgelist, paralogs_wgd, gma_ppi)
    motifs_count <- find_delta(
        edgelist, paralogs_wgd, gma_ppi, count_only = TRUE
    )
    
    # Simulate presence of delta motif for testing purposes
    fake_grn <- data.frame(
        Node1 = c("Node1", "Node1"),
        Node2 = c("Node2", "Node3")
    )
    fake_paralogs <- data.frame(
        Node1 = "Node2",
        Node2 = "Node3"
    )
    fake_ppi <- data.frame(
        Node1 = "Node2",
        Node2 = "Node3"
    )
    motif_fake <- find_delta(fake_grn, fake_paralogs, fake_ppi)
    
    expect_true(is.null(motifs))
    expect_true(is.numeric(motifs_count))
    expect_equal(length(motifs_count), 1)
    expect_equal(class(motif_fake), "character")
    expect_equal(length(motif_fake), 1)
    
})

test_that("find_ppi_v() returns motifs as a character vector", {
    edgelist <- gma_ppi
    motifs <- find_v(edgelist, paralogs_wgd)
    expect_equal(class(motifs), "character")
})

test_that("find_v() returns motifs as a character vector", {
    edgelist <- gma_grn[1:4000, 1:2] # reducing for test purposes
    motifs <- find_v(edgelist, paralogs_wgd)
    expect_true(class(motifs) %in% c("character", "NULL"))
})

test_that("find_bifan() returns motifs as an edge list", {
    paralogs <- rbind(
        paralogs_wgd,
        data.frame(duplicate1 = "Glyma.01G177200", 
                   duplicate2 = "Glyma.08G116700")
    )
    
    motifs <- find_bifan(gma_grn[600:1400, ], paralogs)
    motifs_count <- find_bifan(gma_grn[600:1400, ], paralogs, count_only = TRUE)
    
    expect_equal(class(motifs), "character")
    expect_true(is.numeric(motifs_count))
    expect_equal(length(motifs_count), 1)
})

