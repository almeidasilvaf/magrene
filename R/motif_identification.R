
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


#' Check if the proteins encoded by two genes interact
#' 
#' @param gene1 Character of gene ID for gene 1.
#' @param gene2 Character of gene ID for gene 2.
#' @param edgelist_ppi A 2-column data frame with IDs of genes that encode
#' each protein in the interacting pair.
#' 
#' @return Logical indicating whether the proteins encoded 
#' by \strong{gene1} and \strong{gene2} interact in a PPI network.
#' 
#' @export
#' @rdname are_interacting
#' @examples 
#' data(gma_ppi)
#' gene1 <- "Glyma.19G213200"
#' gene2 <- "Glyma.01G004300"
#' are_paralogs(gene1, gene2, gma_ppi)
are_interacting <- function(gene1 = NULL, gene2 = NULL, edgelist_ppi = NULL) {
    
    names(edgelist_ppi) <- c("gene1", "gene2")
    partners_gene1 <- edgelist_ppi[edgelist_ppi$gene1 == gene1 |
                                       edgelist_ppi$gene2 == gene1, ]
    partners_gene1 <- unique(c(as.character(partners_gene1$gene1),
                               as.character(partners_gene1$gene2)))
    status <- FALSE
    if(gene2 %in% partners_gene1) {
        status <- TRUE
    }
    return(status)
}

#' Find lambda motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#'
#' @return A list of data frames, each data frame being an egde list with
#' the interactions in the motif.
#' @importFrom utils combn
#' @export
#' @rdname find_lambda
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' motifs <- find_lambda(edgelist, paralogs)
find_lambda <- function(edgelist = NULL, paralogs = NULL) {
    
    # Create list of TFs and its interacting partners
    names(edgelist) <- c("Node1", "Node2")
    tfs_and_interactions <- split(edgelist, edgelist$Node1)
    len <- vapply(tfs_and_interactions, nrow, numeric(1))
    tfs_and_interactions <- tfs_and_interactions[len >= 2]
    
    # Create a list of all combinations of partners for each TF, then
    # count how many combinations include a paralog pair
    motifs <- lapply(tfs_and_interactions, function(x) {
        partners <- unique(x[, 2])
        comb <- utils::combn(partners, 2, simplify = FALSE)
        paralog_partners <- lapply(comb, function(y) {
            check_paralogs <- are_paralogs(y[1], y[2], paralogs)
            return(check_paralogs)
        })
        lambda.idx <- which(unlist(paralog_partners) == TRUE)
        if(length(lambda.idx) >= 1) {
            edges <- lapply(lambda.idx, function(i) {
                target_nodes <- comb[[i]]
                e <- x[x$Node2 %in% target_nodes, ]
                return(e)
            })
        } else {
            edges <- NULL
        }
        return(edges)
    })
    remove <- vapply(motifs, is.null, logical(1))
    final_motifs <- unlist(motifs[!remove], recursive = FALSE)
    return(final_motifs)
}


#' Find delta motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param edgelist_ppi A 2-column data frame with IDs of genes that encode
#' each protein in the interacting pair.
#' 
#' @return A list of data frames, each data frame being an egde list with
#' the interactions in the motif.
#' @importFrom utils combn
#' @export
#' @rdname find_delta
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' data(gma_ppi)
#' edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' edgelist_ppi <- gma_ppi
#' motifs <- find_delta(edgelist, paralogs, edgelist_ppi)
find_delta <- function(edgelist = NULL, paralogs = NULL,
                       edgelist_ppi = NULL) {
    
    # Create list of TFs and its interacting partners
    names(edgelist) <- c("Node1", "Node2")
    tfs_and_interactions <- split(edgelist, edgelist$Node1)
    len <- vapply(tfs_and_interactions, nrow, numeric(1))
    tfs_and_interactions <- tfs_and_interactions[len >= 2]
    
    # Create a list of all combinations of partners for each TF, then
    # count how many combinations include a paralog pair
    motifs <- lapply(tfs_and_interactions, function(x) {
        partners <- unique(x[, 2])
        comb <- utils::combn(partners, 2, simplify = FALSE)
        paralog_partners <- lapply(comb, function(y) {
            check_paralogs <- are_paralogs(y[1], y[2], paralogs)
            check_ppi <- FALSE
            if(check_paralogs) {
                check_ppi <- are_interacting(y[1], y[2], edgelist_ppi)
            }
            check <- FALSE
            if(check_paralogs == TRUE & check_ppi == TRUE) {
                check <- TRUE
            }
            return(check)
        })
        lambda.idx <- which(unlist(paralog_partners) == TRUE)
        if(length(lambda.idx) >= 1) {
            edges <- lapply(lambda.idx, function(i) {
                target_nodes <- comb[[i]]
                e <- x[x$Node2 %in% target_nodes, ]
                return(e)
            })
        } else {
            edges <- NULL
        }
        return(edges)
    })
    remove <- vapply(motifs, is.null, logical(1))
    final_motifs <- unlist(motifs[!remove], recursive = FALSE)
    return(final_motifs)
}


#' Find V motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#'
#' @return A list of data frames, each data frame being an egde list with
#' the interactions in the motif.
#' @importFrom utils combn
#' @export
#' @rdname find_v
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[2000:4000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' motifs <- find_v(edgelist, paralogs)
find_v <- function(edgelist = NULL, paralogs = NULL) {
    
    # Create list of TFs and its interacting partners
    names(edgelist) <- c("Node1", "Node2")
    targets_and_tfs <- split(edgelist, edgelist$Node2)
    len <- vapply(targets_and_tfs, nrow, numeric(1))
    targets_and_tfs <- targets_and_tfs[len >= 2]
    
    # Create a list of all combinations of partners for each TF, then
    # count how many combinations include a paralog pair
    motifs <- lapply(targets_and_tfs, function(x) {
        tfs <- unique(x$Node1)
        comb <- utils::combn(tfs, 2, simplify = FALSE)
        paralog_tfs <- lapply(comb, function(y) {
            check_paralogs <- are_paralogs(y[1], y[2], paralogs)
            return(check_paralogs)
        })
        v.idx <- which(unlist(paralog_tfs) == TRUE)
        if(length(v.idx) >= 1) {
            edges <- lapply(v.idx, function(i) {
                tf_nodes <- comb[[i]]
                e <- x[x$Node1 %in% tf_nodes, ]
                return(e)
            })
        } else {
            edges <- NULL
        }
        return(edges)
    })
    remove <- vapply(motifs, is.null, logical(1))
    final_motifs <- unlist(motifs[!remove], recursive = FALSE)
    return(final_motifs)
}







