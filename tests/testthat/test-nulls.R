
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
    n <- 2
    nulls <- generate_nulls(edgelist[1:100, ], paralogs, edgelist_ppi, n)
    expect_equal(class(nulls), "list")
    expect_true(class(nulls[[1]]) %in% c("integer", "numeric"))
    expect_true("lambda" %in% names(nulls))
    expect_true("delta" %in% names(nulls))
    expect_true("V" %in% names(nulls))
    expect_true("bifan" %in% names(nulls))
    expect_true("PPI_V" %in% names(nulls))
})

test_that("calculate_Z() returns a numeric vector of Z-scores", {
    null <- rnorm(1000, mean = 5, sd = 1)
    nulls <- list(lambda = null, V = null, delta = null, bifan = null)
    observed <- list(lambda = 7, bifan = 13, delta = 9, V = 5)
    z <- calculate_Z(observed, nulls)
    expect_equal(class(z), "numeric")
    expect_equal(length(z), 4)
})
