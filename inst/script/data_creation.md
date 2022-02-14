
# Data in data/

## gma\_grn.rda

BioProjects:

-   PRJNA79597 (Libault *et al.*, 2010. The Plant Journal)
-   PRJNA208048 (Severin *et al.*, 2010. BMC Plant Biology)

The gene expression matrix (TPM, kallisto) and sample metadata were
downloaded from the Soybean Expression Atlas (Machado *et al.*, 2020).

``` r
load("~/Dropbox/Atlas/atlasv2_tpm.rda")
gmax_bioproj <- c("PRJNA79597", "PRJNA208048")
gmax_exp <- atlas_tpm[, atlas_tpm$BioProject %in% gmax_bioproj]
gmax_exp <- BioNERO::exp_preprocess(
    gmax_exp, 
    min_exp = 5, 
    variance_filter = TRUE, n = 4000
)
rm(atlas_tpm)
```

Transcription factors:

``` r
gmax_tfs <- readr::read_tsv("http://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz", show_col_types = FALSE)[, c("Gene_ID", "Family")]
gmax_tfs <- as.data.frame(gmax_tfs)

reg <- unique(gmax_tfs$Gene_ID)
```

GRN inference:

``` r
library(BioNERO)
set.seed(123)
genie3 <- BioNERO::grn_infer(
    gmax_exp, 
    method = "genie3", 
    regulators = reg
)

gma_grn <- grn_filter(genie3, nsplit = 20)
# 137866
gma_grn <- genie3[1:137866, ]

usethis::use_data(gma_grn, compress = "xz")
```

## gma\_paralogs.rda

Paralogs were downloaded from the Supplementary Data of [Almeida-Silva
*et al.*, 2020](https://doi.org/10.1007/s00425-020-03499-8).

``` r
files <- c(
    TD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/td_kaks.txt",
    PD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/pd_kaks.txt",
    TRD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/trd_kaks.txt",
    WGD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/wgd_kaks.txt",
    DD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/dd_kaks.txt"
)

gma_paralogs <- Reduce(rbind, lapply(seq_along(files), function(x) {
    df <- readr::read_tsv(files[[x]], col_names = TRUE, 
                          show_col_types = FALSE)[-1, 1:2]
    df <- cbind(df, Type = names(files)[x])
    return(df)
}))
names(gma_paralogs) <- c("duplicate1", "duplicate2", "type")
usethis::use_data(gma_paralogs, compress = "xz")
```

## gma\_ppi.rda

Retrieve genes included in the GRN.

``` r
data(gma_grn)
genes <- unique(c(gma_grn$Node1, gma_grn$Node2))
genes <- gsub("Glyma\\.", "GLYMA_", genes)
```

The PPI network and protein ID correspondences were downloaded from
STRING.

``` r
# Create a table of protein to gene correspondence
gma_aliases <- read.csv(
    "~/Downloads/3847.protein.aliases.v11.5.txt.gz", sep = "\t", header = TRUE
)

gma_ids <- gma_aliases[gma_aliases$source == "Ensembl_UniProt_GN", ]
gma_ids <- gma_ids[startsWith(gma_ids$alias, "GLYMA"), 1:2]
gma_ids$gene <- gsub("GLYMA_", "Glyma.", gma_ids$alias)
gma_ids <- gma_ids[, -2]
names(gma_ids) <- c("STRING_ID", "gene")
rm(gma_aliases)
```

Create the final PPI.

``` bash
# Keep only interactions with confidence scores > 0.4
cd ~/Downloads
cat 3847.protein.physical.links.v11.5.txt.gz | zcat | awk '($3 > 400) {print}' > gma_ppi.txt
```

``` r
# Read files and convert protein names to Glyma IDs
gma_ppi <- read.csv(
    "~/Downloads/gma_ppi.txt", header = TRUE, sep = " "
)[, 1:2]
gma_ppi <- merge(gma_ppi, gma_ids, by.x = "protein1", by.y = "STRING_ID")
gma_ppi$protein1 <- NULL
gma_ppi <- merge(gma_ppi, gma_ids, by.x = "protein2", by.y = "STRING_ID")
gma_ppi$protein2 <- NULL
names(gma_ppi) <- c("node1", "node2")


# Filter PPI network to only keep genes in the GRN
genes_grn <- unique(c(gma_grn$Node1, gma_grn$Node2))
gma_ppi <- gma_ppi[gma_ppi$node1 %in% genes_grn, ]
gma_ppi <- gma_ppi[gma_ppi$node2 %in% genes_grn, ]
usethis::use_data(gma_ppi, compress = "xz")
```
