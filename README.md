MultiEWCE
================
[![](https://img.shields.io/badge/devel%20version-0.1.3-black.svg)](https://github.com/neurogenomics/MultiEWCE)
[![R build
status](https://github.com/neurogenomics/MultiEWCE/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/MultiEWCE/actions)
[![](https://img.shields.io/github/last-commit/neurogenomics/MultiEWCE.svg)](https://github.com/neurogenomics/MultiEWCE/commits/master)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/MultiEWCE.svg)](https://github.com/neurogenomics/MultiEWCE)
[![](https://codecov.io/gh/neurogenomics/MultiEWCE/branch/master/graph/badge.svg)](https://codecov.io/gh/neurogenomics/MultiEWCE)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
¶ <h4> ¶ Authors: <i>Robert Gordon-Smith, Brian Schilder, Nathan
Skene</i> ¶ </h4>
<h4> ¶ README updated: <i>Feb-09-2023</i> ¶ </h4>

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

> Nathan G. Skene, Seth G. N. Grant (2016) Identification of cell types
> underlying over 8,000 rare diseases and subtraits, *bioRxiv*;

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
