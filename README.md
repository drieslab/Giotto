
<!-- README.md is generated from README.Rmd. Please edit that file -->
Giotto
======

<!-- badges: start -->
<!-- badges: end -->
The goal of Giotto is to process, analyze and visualize **single-cell spatial transcriptomic** data. Simultaneously this package contains the data that was used in the recent [**seqFISH+**](https://www.nature.com/articles/s41586-019-1049-y) paper and can thus be used to explore or re-analyze this dataset.

Requirements
------------

-   R (&gt;= 3.5.1)
-   Python (&gt;= 3.0)
-   Unix/Linux

Installation
------------

#### R installation

You can install the development version of Giotto with:

``` r
library(devtools)
install_github("RubD/Giotto")
```

You can NOT YET install the released version of Giotto from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Giotto")
```

#### Python tools (optional)

This is necessary to run all available analyses, including Leiden/Louvain clustering and [**HMRF**](http://www.nature.com/articles/nbt.4260) and to build and use the interactive visualization tool. An alternative, but less flexible, R version for Louvain clustering is available.

``` python
conda install GioTools
```

Examples
--------

<img src="./inst/extdata/cortex_svz_location_fields.png" alt="raw image" width="264" />

[<img src="./inst/images/test_cortex.png" alt="Cortex" width="377" />](./inst/examples/mouse_cortex_svz/mouse_cortex_example.md) <img src="./inst/images/test_OB.png" alt="OB" width="377" />

Latest News
-----------

-   First release (beta)
-   ...

FAQ
---

-   Cortex/SVZ and olfactory bulb data availability?
-   ...

References
----------

-   seqFISH+ paper
-   HMRF paper
