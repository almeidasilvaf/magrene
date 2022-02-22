
#' Generate null distributions of motif counts for each motif type
#' 
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param edgelist_ppi A 2-column data frame with IDs of genes that encode
#' each protein in the interacting pair.
#' @param n Number of degree-preserving simulated networks to generate.
#' Default: 1000.
#'
#' @return A list of numeric vectors named `lambda`, `delta`, `V`, and `bifan`, 
#' containing the null distribution of motif counts for each motif type.
#' 
#' @export
#' @rdname generate_nulls
#' @examples 
#' set.seed(123)
#' data(gma_grn)
#' data(gma_paralogs)
#' data(gma_ppi)
#' edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' edgelist_ppi <- gma_ppi
#' n <- 2 # small n for demonstration purposes
#' generate_nulls(edgelist, paralogs, edgelist_ppi, n)
generate_nulls <- function(edgelist = NULL, paralogs = NULL, 
                           edgelist_ppi = NULL, n = 1000) {
    
    nulls <- lapply(seq_len(n), function(x) { # for each iteration x
        # Simulate network by shuffling target genes
        names(edgelist) <- c("Node1", "Node2")
        sim_grn <- edgelist
        sim_grn$Node2 <- sample(sim_grn$Node2, replace = FALSE)

        # Calculate motif frequencies in iteration x and store them in a vector
        n_lambda <- length(find_lambda(sim_grn, paralogs))
        n_delta <- length(find_delta(sim_grn, paralogs, edgelist_ppi))
        n_v <- length(find_v(sim_grn, paralogs))
        n_bifan <- length(find_bifan(sim_grn, paralogs))
        
        n_iteration <- c(n_lambda, n_delta, n_v, n_bifan)
        names(n_iteration) <- c("lambda", "delta", "V", "bifan")
        return(n_iteration)
    })
    
    # Create vectors of null distros for each motif type
    nulls_vector <- unlist(nulls)
    lambda_distro <- nulls_vector[names(nulls_vector) == "lambda"]
    delta_distro <- nulls_vector[names(nulls_vector) == "delta"]
    v_distro <- nulls_vector[names(nulls_vector) == "V"]
    bifan_distro <- nulls_vector[names(nulls_vector) == "bifan"]
    
    # Store results in a list
    null_list <- list(
        lambda = lambda_distro, delta = delta_distro,
        V = v_distro, bifan = bifan_distro
    )
    return(null_list)
}



#' Calculate Z-score for motif frequencies
#'
#' @param observed A list of observed motif frequencies for each motif type.
#' List elements must be named 'lambda', 'bifan', 'V', and 
#' 'delta' (not necessarily in that order).
#' @param nulls A list of null distributions for each motif type 
#' as returned by \code{generate_nulls}.
#' 
#' @return A numeric vector with the Z-score for each motif type.
#' 
#' @importFrom stats sd
#' @export
#' @rdname calculate_Z
#' @examples
#' # Simulating it for test purposes
#' null <- rnorm(1000, mean = 5, sd = 1)
#' nulls <- list(lambda = null, V = null, delta = null, bifan = null)
#' observed <- list(lambda = 7, bifan = 13, delta = 9, V = 5)
#' z <- calculate_Z(observed, nulls)
#' # Check for motif enrichment (Z > 5)
#' z[which(z > 5)]
calculate_Z <- function(observed = NULL, nulls = NULL) {
    
    Z_score <- unlist(lapply(names(nulls), function(x) {
        s <- ( observed[[x]] - mean(nulls[[x]]) ) / stats::sd(nulls[[x]])
        return(s)
    }))
    names(Z_score) <- names(nulls)
    return(Z_score)
}


