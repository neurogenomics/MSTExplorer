MSTExplorer
================
<img src='https://github.com/neurogenomics/MSTExplorer/raw/master/inst/hex/hex.png' title='Hex sticker for MSTExplorer' height='300'><br>
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/devel%20version-1.0.7-black.svg)](https://github.com/neurogenomics/MSTExplorer)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/MSTExplorer.svg)](https://github.com/neurogenomics/MSTExplorer)
[![](https://img.shields.io/github/last-commit/neurogenomics/MSTExplorer.svg)](https://github.com/neurogenomics/MSTExplorer/commits/master)
<br> [![R build
status](https://github.com/neurogenomics/MSTExplorer/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/MSTExplorer/actions)
[![](https://codecov.io/gh/neurogenomics/MSTExplorer/branch/master/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/MSTExplorer)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/MSTExplorer/tree/master' target='_blank'><img src='https://codecov.io/gh/neurogenomics/MSTExplorer/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Brian Schilder, Robert Gordon-Smith, Nathan Skene,
Hiranyamaya Dash</i>  
</h4>
<h4>  
README updated: <i>Mar-09-2025</i>  
</h4>

<!-- To modify Package/Title/Description/Authors fields, edit the DESCRIPTION file -->

## Introduction

Many genes have been associated with diseases Multi-Scale Target
Explorer (`MSTExplorer`) systematically identifies, prioritises, and
visualises cell-type-specific gene therapy targets across the phenome.

Core functionalities include:

**1. Conducting phenotype x cell type genetic association tests at
scale**

- The [Human Phenotype Ontology](https://hpo.jax.org/) (integrated with
  gene annotations from [OMIM](https://omim.org/) /
  [DECIPHER](https://www.deciphergenomics.org/) /
  [ORPHANET](https://www.orpha.net/)) is used as the source of phenotype
  gene signatures. Each gene-phenotype associated is given a continuous
  score that approximates the current strength of evidence for the
  association (using data derived from [GenCC](https://thegencc.org/)).

- Whole-body scRNA-seq atlases from humans (across multiple
  developmental stages) are used as a data-driven source of cell
  type-specific gene markers.  

- The underlying association tests are designed for both speed and
  accuracy using memory-efficient data structures, and a highly
  parallelisable implementation of Generalised Linear Regression (GLM).
  For example, associations for all pairwise combinations of \>11k
  phenotypes x \>200 cell types (\>2,200,000 associations) can be in
  \<30 minutes on a Macbook laptop with 10 CPU cores).

**2. Inferring multi-scale causal graphs of disease**

`MSTExplorer` allows users to easily infer and construct multi-scale
causal graphs of Diseases (blue nodes) -\> Phenotypes (purple nodes) -\>
Cell types (orange nodes) -\> Genes (yellow nodes).

<figure>
<img
src="https://github.com/neurogenomics/rare_disease_celltyping/blob/299abe0ccd00644bc2f05a1389704fe196a3e868/manuscript/_manuscript/img/fig-therapy-examples-supp/lethal_skeletal_dysplasia.png?raw=true"
height="400"
alt="Example multi-scale network focused on lethal skeletal dysplasia, a phenotype of multiple diseases" />
<figcaption aria-hidden="true"><em>Example multi-scale network focused
on lethal skeletal dysplasia, a phenotype of multiple
diseases</em></figcaption>
</figure>

[See here for more example
networks.](https://github.com/neurogenomics/rare_disease_celltyping/tree/299abe0ccd00644bc2f05a1389704fe196a3e868/manuscript/_manuscript/img/fig-therapy-examples-supp).

**3. Prioritising cell-type-specific gene therapy targets**

`MSTExplorer` also provides a comprehensive and customisable pipeline
that can be run via a single function (`prioritise_targets()`) to
produce the most promising cell-type-specific gene therapy targets
across the phenome.

## Installation

Within R:

``` r
if(!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("neurogenomics/MSTExplorer")
library(MSTExplorer)
```

## Documentation

#### [Website](https://neurogenomics.github.io/MSTExplorer)

#### [Get started](https://neurogenomics.github.io/MSTExplorer/articles/MSTExplorer)

#### [Docker/Singularity Container](https://neurogenomics.github.io/MSTExplorer/articles/docker.html)

## Citation

If you use `MSTExplorer`, please cite:

<!-- Modify this by editing the file: inst/CITATION  -->

> Kitty B. Murphy, Robert Gordon-Smith, Jai Chapman, Momoko Otani, Brian
> M. Schilder, Nathan G. Skene (2023) Identification of cell
> type-specific gene targets underlying thousands of rare diseases and
> subtraits. medRxiv, <https://doi.org/10.1101/2023.02.13.23285820>

## Contact

### [Neurogenomics Lab](https://www.neurogenomics.co.uk)

UK Dementia Research Institute  
Department of Brain Sciences  
Faculty of Medicine  
Imperial College London  
[GitHub](https://github.com/neurogenomics)

## Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.3.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        jsonlite_1.9.1      renv_1.1.2         
    ##  [4] dplyr_1.1.4         compiler_4.4.2      BiocManager_1.30.25
    ##  [7] tidyselect_1.2.1    rvcheck_0.2.1       scales_1.3.0       
    ## [10] yaml_2.3.10         fastmap_1.2.0       here_1.0.1         
    ## [13] ggplot2_3.5.1       R6_2.6.1            generics_0.1.3     
    ## [16] knitr_1.49          yulab.utils_0.2.0   tibble_3.2.1       
    ## [19] desc_1.4.3          dlstats_0.1.7       munsell_0.5.1      
    ## [22] rprojroot_2.0.4     pillar_1.10.1       RColorBrewer_1.1-3 
    ## [25] rlang_1.1.5         badger_0.2.4        xfun_0.51          
    ## [28] fs_1.6.5            cli_3.6.4           magrittr_2.0.3     
    ## [31] rworkflows_1.0.6    digest_0.6.37       grid_4.4.2         
    ## [34] rstudioapi_0.17.1   lifecycle_1.0.4     vctrs_0.6.5        
    ## [37] evaluate_1.0.3      glue_1.8.0          data.table_1.17.0  
    ## [40] colorspace_2.1-1    rmarkdown_2.29      tools_4.4.2        
    ## [43] pkgconfig_2.0.3     htmltools_0.5.8.1

</details>

<hr>
