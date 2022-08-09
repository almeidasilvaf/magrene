
<!-- README.md is generated from README.Rmd. Please edit that file -->

# magrene <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/almeidasilvaf/magrene/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/almeidasilvaf/magrene/actions)
<!-- badges: end -->

The goal of `magrene` is to identify and analyze graph motifs containing
duplicated genes in gene regulatory networks (GRNs). Possible motifs
include V, PPI V, lambda, delta, and bifans.

<img src="vignettes/motifs_vignette.png" title="Network motifs and functions to identify them. Shaded boxes indicate paralogs. Regulators and targets are indicated in purple and green, respectively. Arrows indicate directed regulatory interactions, while dashed lines indicate protein-protein interaction." alt="Network motifs and functions to identify them. Shaded boxes indicate paralogs. Regulators and targets are indicated in purple and green, respectively. Arrows indicate directed regulatory interactions, while dashed lines indicate protein-protein interaction." width="100%" />

Motif frequencies can be analyzed in the context of gene duplications to
explore the impact of small-scale and whole-genome duplications on gene
regulatory networks.

Additionally, GRNs can be tested for motif enrichment by comparing motif
frequencies to a null distribution generated from degree-preserving
simulated GRNs.

Finally, users can calculate the interaction similarity between gene
pairs based on the Sorensen-Dice similarity index.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `magrene` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("magrene")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/magrene) with:

``` r
BiocManager::install("almeidasilvaf/magrene")
```

## Citation

Below is the citation output from using `citation('magrene')` in R.
Please run this yourself to check for any updates on how to cite
**magrene**.

``` r
print(citation('magrene'), bibtex = TRUE)
#> 
#> To cite package 'magrene' in publications use:
#> 
#>   Almeida-Silva F, Van de Peer Y (2022). _magrene: Motif Analysis In
#>   Gene Regulatory Networks_. R package version 0.99.0,
#>   <https://github.com/almeidasilvaf/magrene>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {magrene: Motif Analysis In Gene Regulatory Networks},
#>     author = {Fabrício Almeida-Silva and Yves {Van de Peer}},
#>     year = {2022},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/almeidasilvaf/magrene},
#>   }
```

Please note that the `magrene` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `magrene` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.15/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://almeidasilvaf.github.io/magrene)
    is automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.15/biocthis)*.
