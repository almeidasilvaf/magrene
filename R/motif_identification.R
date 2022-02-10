
#' Check if two genes are paralogs
#' 
#' @param gene1 Character of gene ID for gene 1.
#' @param gene2 Character of gene ID for gene 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' 
#' @return Logical indicating whether \strong{gene1} and \strong{gene2} are
#' paralogs.
#' 
#' @export
#' @rdname are_paralogs
#' @examples 
#' data(gma_paralogs)
#' paralogs <- gma_paralogs[, 1:2]
#' gene1 <- "Glyma.01G011500"
#' gene2 <- "Glyma.01G011600"
#' are_paralogs(gene1, gene2, paralogs)
are_paralogs <- function(gene1 = NULL, gene2 = NULL, paralogs = NULL) {
    
    names(paralogs) <- c("duplicate1", "duplicate2")
    paralogs_gene1 <- paralogs[paralogs$duplicate1 == gene1 |
                                   paralogs$duplicate2 == gene1, ]
    paralogs_gene1 <- unique(c(as.character(paralogs_gene1$duplicate1),
                               as.character(paralogs_gene1$duplicate2)))
    status <- FALSE
    if(gene2 %in% paralogs_gene1) {
        status <- TRUE
    }
    return(status)
}


#' Count the frequency of lambda motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#'
#' @return Numeric scalar with the frequency of lambda motifs in the network.
#' @importFrom utils combn
#' @export
#' @rdname count_lambda
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[1:500, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' count <- count_lambda(edgelist, paralogs)
count_lambda <- function(edgelist = NULL, paralogs = NULL) {
    
    # Create list of TFs and its interacting partners
    tfs_and_interactions <- split(edgelist, edgelist[,1])
    len <- vapply(tfs_and_interactions, nrow, numeric(1))
    tfs_and_interactions <- tfs_and_interactions[len >= 2]
    
    # Create a list of all combinations of partners for each TF
    count <- lapply(tfs_and_interactions, function(x) {
        partners <- unique(x[, 2])
        comb <- utils::combn(partners, 2, simplify = FALSE)
        paralog_partners <- lapply(comb, function(y) {
            check <- are_paralogs(y[1], y[2], paralogs)
            return(check)
        })
        lambda_count <- sum(unlist(paralog_partners))
        return(lambda_count)
    })
    final_count <- sum(unlist(count))
    return(final_count)
}


