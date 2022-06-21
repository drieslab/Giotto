# `readExprMatrix`

Read expression matrix


## Description

Function to read an expression matrix into a sparse matrix.


## Usage

```r
readExprMatrix(path, cores = NA, transpose = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`path`     |     path to the expression matrix
`cores`     |     number of cores to use
`transpose`     |     transpose matrix


## Details

The expression matrix needs to have both unique column names and row names


## Value

sparse matrix


