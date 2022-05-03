
#' Calculate Sorensen-Dice similarity between paralogous gene pairs
#'
#' @param edgelist A 2-column data frame with regulators in column 1 and
#' targets in column 2.
#' @param paralogs A 2-column data frame with gene IDs for each paralog
#' in the paralog pair.
#' 
#' @return A data frame containing the paralogous gene pairs and their 
#' Sorensen-Dice similarity scores.
#' @export
#' @rdname sd_similarity
#' @examples 
#' data(gma_ppi)
#' data(gma_paralogs)
#' edgelist <- gma_ppi
#' paralogs <- gma_paralogs
#' sim <- sd_similarity(edgelist, paralogs)
sd_similarity <- function(edgelist = NULL, paralogs = NULL) {
    
    # Create a list of partners for each node
    pl1 <- split(edgelist, edgelist[, 1])
    pl1 <- lapply(pl1, function(x) return(unique(x[,2])))
    
    pl2 <- split(edgelist, edgelist[, 2])
    pl2 <- lapply(pl2, function(x) return(unique(x[,1])))
    
    pl <- c(pl1, pl2)
    partners <- tapply(pl, names(pl), function(x) unlist(x, FALSE, FALSE))
    partners <- lapply(partners, unique)
    
    fpara <- paralogs[paralogs[,1] %in% names(partners) &
                          paralogs[,2] %in% names(partners), ]
    
    # Calculate Sorensen-Dice similarity for pairs
    sim <- unlist(lapply(seq_len(nrow(fpara)), function(x) {
        n1 <- fpara[x, 1]
        n2 <- fpara[x, 2]
        p1 <- partners[[n1]]
        p2 <- partners[[n2]]
        inter <- length(intersect(p1, p2))
        
        sd <- 2 * (inter) / (length(p1) + length(p2))
        return(round(sd, 2))
    }))
    fpara$sorensen_dice <- sim
    return(fpara)
}