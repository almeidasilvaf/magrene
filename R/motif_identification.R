

#' Find lambda motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam()
#' @param count_only Logical indicating whether the function should return
#' only motif counts as a numeric scalar. If FALSE, it will return
#' a character vector of motifs. Default: FALSE.
#'
#' @return A character vector with lambda motifs represented 
#' in the format \strong{target1<-regulator->target2}.
#' @importFrom BiocParallel bplapply SerialParam
#' @export
#' @rdname find_lambda
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' motifs <- find_lambda(edgelist, paralogs)
find_lambda <- function(edgelist = NULL, paralogs = NULL,
                        bp_param = BiocParallel::SerialParam(),
                        count_only = FALSE) {
    
    names(edgelist) <- c("Node1", "Node2")
    names(paralogs) <- c("Dup1", "Dup2")
    
    # Removing edges where target is not in paralogs df
    all_paralogs <- unique(c(paralogs$Dup1, paralogs$Dup2))
    edgelist <- edgelist[edgelist$Node2 %in% all_paralogs, ]
    
    # Filter paralogs df to keep only rows where both duplicates are targets
    paralogs <- paralogs[paralogs$Dup1 %in% edgelist$Node2 &
                             paralogs$Dup2 %in% edgelist$Node2, ]
    motifs <- NULL
    if(nrow(paralogs) > 0) {
        # Create list of all TFs that regulate each target
        reg <- split(edgelist, edgelist$Node2)
        
        # For each paralogous pair, get TFs that regulate both
        motifs <- unlist(bplapply(seq_len(nrow(paralogs)), function(x) {
            t1 <- paralogs[x, 1]
            t2 <- paralogs[x, 2]
            tfs <- intersect(reg[[t1]]$Node1, reg[[t2]]$Node1)
            if(count_only) {
                m <- length(tfs)
            } else {
                m <- NULL
                if(length(tfs) > 0) {
                    m <- unlist(lapply(tfs, function(y) {
                        return(paste0(t1, "<-", y, "->", t2))
                    }))
                }
            }
            return(m)
        }, BPPARAM = bp_param))
    }
    
    if(count_only) {
        motifs <- sum(motifs)
    }
    return(motifs)
}


#' Find delta motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2. It can be ignored if you give lambda motifs
#' to parameter \strong{lambda_vec} (recommended).
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair. It can be ignored if you give lambda motifs
#' to parameter \strong{lambda_vec} (recommended).
#' @param edgelist_ppi A 2-column data frame with IDs of genes that encode
#' each protein in the interacting pair.
#' @param lambda_vec A character of lambda motifs as returned 
#' by \code{find_lambda()}. If this is NULL, this function will find lambda
#' motifs from \strong{edgelist} and \strong{paralogs} first. Passing 
#' previously identified lambda motifs will make this function much faster.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam().
#' @param count_only Logical indicating whether the function should return
#' only motif counts as a numeric scalar. If FALSE, it will return
#' a character vector of motifs. Default: FALSE.
#' 
#' @return A character vector with lambda motifs represented 
#' in the format \strong{target1<-regulator->target2}.
#' @export
#' @rdname find_delta
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' data(gma_ppi)
#' edgelist <- gma_grn[500:1000, 1:2] # reducing for test purposes
#' edgelist <- gma_grn[1:10000, 1:2]
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' edgelist_ppi <- gma_ppi
#' lambda_vec <- find_lambda(edgelist, paralogs)
#' motifs <- find_delta(edgelist_ppi = edgelist_ppi, lambda_vec = lambda_vec)
find_delta <- function(edgelist = NULL, paralogs = NULL,
                       edgelist_ppi = NULL, lambda_vec = NULL,
                       bp_param = BiocParallel::SerialParam(),
                       count_only = FALSE) {
    
    ivec <- paste0(edgelist_ppi[, 1], edgelist_ppi[, 2]) # PPI vector
    
    # Look for lambda motifs for which targets interact (delta)
    lambda <- lambda_vec
    if(is.null(lambda_vec)) {
        lambda <- find_lambda(edgelist, paralogs)
    }
    motifs <- unlist(bplapply(lambda, function(x) {
        t1 <- gsub("<-.*", "", x)
        t2 <- gsub(".*->", "", x)
        
        pair1 <- paste0(t1, t2)
        pair2 <- paste0(t2, t1)
        check1 <- pair1 %in% ivec
        check2 <- pair2 %in% ivec
        check <- check1 + check2
        motif <- NULL
        if(check != 0) {
            motif <- x
        }
    }, BPPARAM = bp_param))
    
    if(count_only) {
        motifs <- length(motifs)
    }
    return(motifs)
}


#' Find V motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam().
#' @param count_only Logical indicating whether the function should return
#' only motif counts as a numeric scalar. If FALSE, it will return
#' a character vector of motifs. Default: FALSE.
#'
#' @return A character vector with V motifs represented 
#' in the format \strong{regulator1->target<-regulator2}.
#' 
#' @importFrom BiocParallel bplapply SerialParam
#' @export
#' @rdname find_v
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[2000:4000, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' motifs <- find_v(edgelist, paralogs)
find_v <- function(edgelist = NULL, paralogs = NULL,
                   bp_param = BiocParallel::SerialParam(),
                   count_only = FALSE) {
    
    names(edgelist) <- c("Node1", "Node2")
    names(paralogs) <- c("Dup1", "Dup2")
    
    # Removing edges where regulator is not in paralogs df
    all_paralogs <- unique(c(paralogs$Dup1, paralogs$Dup2))
    edgelist <- edgelist[edgelist$Node1 %in% all_paralogs, ]
    
    # Filter paralogs df to keep only rows where both duplicates are TFs
    paralogs <- paralogs[paralogs$Dup1 %in% edgelist$Node1 &
                             paralogs$Dup2 %in% edgelist$Node1, ]
    motifs <- NULL
    if(nrow(paralogs) > 0) {
        # Create list of all targets for each regulator
        tar <- split(edgelist, edgelist$Node1)
        
        # For each paralogous pair, get TFs that regulate both
        motifs <- unlist(bplapply(seq_len(nrow(paralogs)), function(x) {
            t1 <- paralogs[x, 1]
            t2 <- paralogs[x, 2]
            tars <- intersect(tar[[t1]]$Node2, tar[[t2]]$Node2)
            if(count_only) {
                m <- length(tars)
            } else {
                m <- NULL
                if(length(tars) > 0) {
                    m <- unlist(lapply(tars, function(y) {
                        return(paste0(t1, "<-", y, "->", t2))
                    }))
                }
            }
            return(m)
        }, BPPARAM = bp_param))
    }
    if(count_only) {
        motifs <- sum(motifs)
    }
    return(motifs)
}


#' Find V motifs in protein-protein interactions
#'
#' @param edgelist A 2-column data frame with protein 1 in column 1 and
#' protein 2 in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam().
#' @param count_only Logical indicating whether the function should return
#' only motif counts as a numeric scalar. If FALSE, it will return
#' a character vector of motifs. Default: FALSE.
#'
#' @return A character vector with V motifs represented 
#' in the format \strong{paralog1-partner-paralog2}.
#' 
#' @details This function aims to find the number of paralogous gene
#' pairs that share an interaction partner.
#' 
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom utils combn
#' @export
#' @rdname find_ppi_v
#' @examples 
#' data(gma_ppi)
#' data(gma_paralogs)
#' edgelist <- gma_ppi
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' motifs <- find_ppi_v(edgelist, paralogs)
find_ppi_v <- function(edgelist = NULL, paralogs = NULL,
                       bp_param = BiocParallel::SerialParam(),
                       count_only = FALSE) {
    
    names(edgelist) <- c("Node1", "Node2")
    names(paralogs) <- c("Dup1", "Dup2")
    
    # Removing edges where target is not in paralogs df
    all_paralogs <- unique(c(paralogs$Dup1, paralogs$Dup2))
    edgelist <- edgelist[edgelist$Node1 %in% all_paralogs, ]
     
    pvec <- paste0(paralogs$Dup1, paralogs$Dup2) # paralogs vector
    
    # Create list of proteins and their interacting partners
    targets_and_p <- split(edgelist, edgelist$Node2)
    len <- vapply(targets_and_p, nrow, numeric(1))
    targets_and_p <- targets_and_p[len >= 2]
    
    motifs <- NULL
    if(length(targets_and_p) > 0) {
        motifs <- unlist(bplapply(targets_and_p, function(x) {
            p <- unique(x[, 1])
            if(length(p) > 1) {
                comb <- utils::combn(p, 2, simplify = FALSE)
                comb1 <- unlist(lapply(comb, function(x) paste0(x[1], x[2])))
                comb2 <- unlist(lapply(comb, function(x) paste0(x[2], x[1])))
                para1 <- comb1 %in% pvec
                para2 <- comb2 %in% pvec
                paralog_partners <- para1 + para2
                v.idx <- which(paralog_partners >= 1)
                if(length(v.idx) >= 1) {
                    edges <- unlist(lapply(v.idx, function(i) {
                        ps <- comb[[i]]
                        e <- paste0(ps[1], "-", x[1,2], "-", ps[2])
                        return(e)
                    }))
                } else {
                    edges <- NULL
                }
            } else {
                edges <- NULL
            }
            return(edges)
        }, BPPARAM = bp_param))
    }
    names(motifs) <- NULL
    if(count_only) { motifs <- length(motifs) }
    return(motifs)
}


#' Find bifan motifs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2. It can be ignored if you give lambda motifs
#' to parameter \strong{lambda_vec} (recommended).
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' @param lambda_vec A character of lambda motifs as returned 
#' by \code{find_lambda()}. If this is NULL, this function will find lambda
#' motifs from \strong{edgelist} and \strong{paralogs} first. Passing 
#' previously identified lambda motifs will make this function much faster.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam().
#' @param count_only Logical indicating whether the function should return
#' only motif counts as a numeric scalar. If FALSE, it will return
#' a character vector of motifs. Default: FALSE.
#'
#' @return A character vector with bifan motifs represented 
#' in the format \strong{regulator1, regulator2->target1, target2}.
#' 
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom utils combn
#' @export
#' @rdname find_bifan
#' @examples 
#' data(gma_grn)
#' data(gma_paralogs)
#' edgelist <- gma_grn[600:1400, 1:2] # reducing for test purposes
#' paralogs <- gma_paralogs[gma_paralogs$type == "WGD", 1:2]
#' paralogs <- rbind(
#' paralogs,
#' data.frame(duplicate1 = "Glyma.01G177200", 
#'            duplicate2 = "Glyma.08G116700")
#' )
#' bifan <- find_bifan(edgelist, paralogs)
find_bifan <- function(edgelist = NULL, paralogs = NULL,
                       lambda_vec = NULL,
                       bp_param = BiocParallel::SerialParam(),
                       count_only = FALSE) {
    
    pvec <- paste0(paralogs[,1], paralogs[,2]) # paralogs vector
    
    # Find lambda motifs and compare them pairwise to look for shared neighbors
    lambdas <- lambda_vec
    if(is.null(lambda_vec)) {
        lambdas <- find_lambda(edgelist, paralogs)
    }
    
    tfnames <- gsub(".*<-", "", lambdas)
    tfnames <- gsub("->.*", "", tfnames)
    final_motifs <- NULL
    if(length(tfnames) > 1) {
        comb <- utils::combn(seq_along(lambdas), 2, simplify = FALSE)
        #comb <- utils::combn(tfnames, 2, simplify = FALSE)
        final_motifs <- unlist(bplapply(comb, function(x) {
            motif <- NULL
            tf1 <- tfnames[x[1]]
            tf2 <- tfnames[x[2]]
            # Check 1) - Are TFs in the lambda motifs paralogs?
            check1 <- paste0(tf1, tf2) %in% pvec
            check2 <- paste0(tf2, tf1) %in% pvec
            check_paralogs <- check1 + check2
            
            # Check 2) Do paralogous TFs regulate the same targets?
            if(check_paralogs != 0) {
                ptf1 <- lambdas[x[1]]
                ptf1 <- c(
                    gsub("<-.*", "", ptf1), gsub(".*->", "", ptf1)
                )
                
                ptf2 <- lambdas[x[2]]
                ptf2 <- c(
                    gsub("<-.*", "", ptf2), gsub(".*->", "", ptf2)
                )
                check_shared <- sum(ptf1 %in% ptf2) == 2
                if(check_shared) {
                    motif <- paste0(tf1, ",", tf2, "->", ptf1[1], ",", ptf1[2])
                }
            }
            return(motif)
        }, BPPARAM = bp_param))
    }
    if(count_only) { final_motifs <- length(final_motifs) }
    return(final_motifs)
}




