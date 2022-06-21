# `read_expression_data`

Read expression data


## Description

Read expression data


## Usage

```r
read_expression_data(
  expr_list = NULL,
  sparse = TRUE,
  cores = NA,
  default_feat_type = NULL,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`expr_list`     |     (nested) list with expression data
`sparse`     |     read matrix data in a sparse manner
`cores`     |     number of cores to use
`default_feat_type`     |     default feature type if nothing is provided
`verbose`     |     be verbose


## Details

mylistA = list('a' = matrix(1:5), 'b' = matrix(1:5))
 depth(mylistA)
 
 mylistB = list(A = list('a' = matrix(1:5), 'b' = matrix(1:5)),
 B = list('c' = matrix(1:5),'d' = matrix(1:5)))
 depth(mylistB)
 
 mylistC = list('RNA' = list('RAW' = list('cell' = matrix(1:5), 'nucleus' = matrix(6:10)),
 'NORM' = list('cell' = matrix(11:15),'nucleus' = matrix(20:25))),
 'PROT' = list('RAW' = list('cell' = matrix(16:20))))
 depth(mylistC)
 
 mymatD = matrix(data = 1:4)


