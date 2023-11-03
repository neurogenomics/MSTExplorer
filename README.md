MultiEWCE
================
NULL [![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/devel%20version-0.1.7-black.svg)](https://github.com/neurogenomics/MultiEWCE)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/MultiEWCE.svg)](https://github.com/neurogenomics/MultiEWCE)
[![](https://img.shields.io/github/last-commit/neurogenomics/MultiEWCE.svg)](https://github.com/neurogenomics/MultiEWCE/commits/master)
<br> [![R build
status](https://github.com/neurogenomics/MultiEWCE/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/MultiEWCE/actions)
[![](https://codecov.io/gh/neurogenomics/MultiEWCE/branch/master/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/MultiEWCE)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/MultiEWCE/tree/master' target='_blank'><img src='https://codecov.io/gh/neurogenomics/MultiEWCE/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Robert Gordon-Smith, Brian Schilder, Nathan Skene</i>  
</h4>
<h4>  
README updated: <i>Nov-03-2023</i>  
</h4>

<!-- To modify Package/Title/Description/Authors fields, edit the DESCRIPTION file -->

## `MultiEWCE`: EWCE for Multiple Gene Lists

### This package allows you to run EWCE on multiple gene lists in parallel.

It has functions for processing the data and running the analyss.
Finally you can merge the results in to one large dataframe. Each
individual gene list analysis will also be saved in the results
directory to allow you to pause the analysis mid way and not lose
anything.

If you use `MultiEWCE`, please cite:

<!-- Modify this by editing the file: inst/CITATION  -->

> Kitty B. Murphy, Robert Gordon-Smith, Jai Chapman, Momoko Otani, Brian
> M. Schilder, Nathan G. Skene (2023) Identification of cell
> type-specific gene targets underlying thousands of rare diseases and
> subtraits. medRxiv, <https://doi.org/10.1101/2023.02.13.23285820>

## Installation

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("neurogenomics/MultiEWCE")
library(MultiEWCE)
```

## Documentation

### [Website](https://neurogenomics.github.io/MultiEWCE)

### [Get started](https://neurogenomics.github.io/MultiEWCE/articles/MultiEWCE)

<br>
