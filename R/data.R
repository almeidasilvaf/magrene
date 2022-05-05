
#' Sample soybean GRN
#' 
#' The GRN was inferred with {BioNERO} using expression data from 
#' Libault et al., 2010, and Severin et al., 2010.
#'
#' @name gma_grn
#' @format A 3-column data frame with node1, node2, and edge weight.
#'
#' @references 
#' Severin, A. J., Woody, J. L., Bolon, Y. T., Joseph, B., Diers, B. W., 
#' Farmer, A. D., ... & Shoemaker, R. C. (2010). 
#' RNA-Seq Atlas of Glycine max: a guide to the soybean transcriptome. 
#' BMC plant biology, 10(1), 1-16.
#' @references
#' Libault, M., Farmer, A., Joshi, T., Takahashi, K., Langley, R. J., 
#' Franklin, L. D., ... & Stacey, G. (2010). 
#' An integrated transcriptome atlas of the crop model Glycine max, 
#' and its use in comparative analyses in plants. 
#' The Plant Journal, 63(1), 86-99.
#' @examples 
#' data(gma_grn)
#' @usage data(gma_grn)
"gma_grn"


#' Soybean (Glycine max) duplicated genes
#' 
#' The repertoire of soybean paralogs was retrieved from 
#' Almeida-Silva et al., 2020.
#'
#' @name gma_paralogs
#' @format A 3-column data frame with duplicate 1, duplicate 2, and 
#' duplication type
#'
#' @references 
#' Almeida-Silva, F., Moharana, K. C., Machado, F. B., & 
#' Venancio, T. M. (2020). Exploring the complexity of soybean (Glycine max) 
#' transcriptional regulation using global gene co-expression networks. 
#' Planta, 252(6), 1-12.
#' @examples 
#' data(gma_paralogs)
#' @usage data(gma_paralogs)
"gma_paralogs"


#' Sample soybean PPI network
#' 
#' PPI were retrieved from the STRING database and filtered to keep only
#' medium confidence edges and nodes in the GRN.
#'
#' @name gma_ppi
#' @format A 2-column data frame with node1 and node2.
#'
#' @examples 
#' data(gma_ppi)
#' @usage data(gma_ppi)
"gma_ppi"


#' Null distribution of motif frequencies for vignette data set
#' 
#' Data were filtered exactly as demonstrated in the vignette.
#'
#' @name nulls
#' @format A list of numeric vectors with the motif frequencies in each
#' simulated network.
#'
#' @examples 
#' data(nulls)
#' @usage data(nulls)
"nulls"

