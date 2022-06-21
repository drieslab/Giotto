# `runPCA_factominer`

runPCA_factominer


## Description

performs PCA based on the factominer package


## Usage

```r
runPCA_factominer(
  x,
  ncp = 100,
  scale = TRUE,
  rev = FALSE,
  set_seed = TRUE,
  seed_number = 1234,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     matrix or object that can be converted to matrix
`ncp`     |     number of principal components to calculate
`scale`     |     scale features
`rev`     |     reverse PCA
`set_seed`     |     use of seed
`seed_number`     |     seed number to use


## Value

list of eigenvalues, loadings and pca coordinates


