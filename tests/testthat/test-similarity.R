
#----Setup----------------------------------------------------------------------
data(gma_ppi)
data(gma_paralogs)
edgelist <- gma_ppi
paralogs <- gma_paralogs

#----Start tests----------------------------------------------------------------
test_that("sd_similarity() returns a data frame with Sorensen-Dice sims", {
    sim <- sd_similarity(edgelist, paralogs)
    expect_equal(class(sim), "data.frame")
    expect_true("sorensen_dice" %in% names(sim))
    expect_true(ncol(sim) >= 3)
})
