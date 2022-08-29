
#----Load data------------------------------------------------------------------
data(gma_paralogs)
data(gma_grn)
data(gma_ppi)

set.seed(123)
edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
edgelist_ppi <- gma_ppi


#----Start tests----------------------------------------------------------------
test_that("generate_nulls() return a list of numeric vectors", {
    n <- 1
    nulls <- generate_nulls(edgelist[1:100, ], paralogs, edgelist_ppi, n)

    # Simulated network with lambda motifs for testing purposes
    fake_grn <- data.frame(
        Node1 = rep(c("Glyma.01G017500", "Glyma.01G017500"), 50),
        Node2 = rep(c("Glyma.01G017200", "Glyma.01G017300"), 50)
    )
    fake_ppi <- data.frame(
        Node1 = rep("Glyma.01G017200", 20),
        Node2 = rep("Glyma.01G017300", 20)
    )
    fake_paralogs <- fake_ppi
    nulls2 <- generate_nulls(fake_grn, fake_paralogs, fake_ppi, n)
    
    # Tests for subset of GRN
    expect_equal(class(nulls), "list")
    expect_true(class(nulls[[1]]) %in% c("integer", "numeric"))
    expect_true("lambda" %in% names(nulls))
    expect_true("delta" %in% names(nulls))
    expect_true("V" %in% names(nulls))
    expect_true("bifan" %in% names(nulls))
    expect_true("PPI_V" %in% names(nulls))
    
    # Tests for fake GRN with lambdas
    expect_equal(class(nulls2), "list")
    expect_true(class(nulls2[[1]]) %in% c("integer", "numeric"))
    expect_true("lambda" %in% names(nulls2))
    expect_true("delta" %in% names(nulls2))
    expect_true("V" %in% names(nulls2))
    expect_true("bifan" %in% names(nulls2))
    expect_true("PPI_V" %in% names(nulls2))
})

test_that("calculate_Z() returns a numeric vector of Z-scores", {
    null <- rnorm(1000, mean = 5, sd = 1)
    nulls <- list(lambda = null, V = null, delta = null, bifan = null)
    observed <- list(lambda = 7, bifan = 13, delta = 9, V = 5)
    z <- calculate_Z(observed, nulls)
    expect_equal(class(z), "numeric")
    expect_equal(length(z), 4)
})
