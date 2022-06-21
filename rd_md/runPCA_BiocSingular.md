# `runPCA_BiocSingular`

runPCA_BiocSingular


## Description

Performs PCA based on the biocSingular package


## Usage

```r
runPCA_BiocSingular(
  x,
  ncp = 100,
  center = TRUE,
  scale = TRUE,
  rev = FALSE,
  set_seed = TRUE,
  seed_number = 1234,
  BSPARAM = c("irlba", "exact", "random"),
  BSParameters = list(NA),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     matrix or object that can be converted to matrix
`ncp`     |     number of principal components to calculate
`center`     |     center the matrix before pca
`scale`     |     scale features
`rev`     |     reverse PCA
`set_seed`     |     use of seed
`seed_number`     |     seed number to use
`BSPARAM`     |     method to use
`BSParameters`     |     additonal parameters for method


## Value

list of eigenvalues, loadings and pca coordinates


