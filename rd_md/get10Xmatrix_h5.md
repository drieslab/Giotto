# `get10Xmatrix_h5`

get10Xmatrix_h5


## Description

This function creates an expression matrix from a 10X h5 file path


## Usage

```r
get10Xmatrix_h5(path_to_data, gene_ids = c("symbols", "ensembl"))
```


## Arguments

Argument      |Description
------------- |----------------
`path_to_data`     |     path to the 10X .h5 file
`gene_ids`     |     use gene symbols (default) or ensembl ids for the gene expression matrix


## Details

If the .h5 10x file has multiple modalities (e.g. RNA and protein),
 multiple matrices will be returned


## Value

(list of) sparse expression matrix from 10X


