# `mini_giotto_single_cell`

mini Giotto object for spatial single-cell resolution data


## Description

Mini Giotto object created from the seqFISH+ data.


## Format

An object of class `"giotto"` ; see [`createGiottoObject`](#creategiottoobject) .


## Usage

```r
data(mini_giotto_single_cell)
```


## References

Eng et al. (2019) Nature
 ( [PubMed](https://www.nature.com/articles/s41586-019-1049-y) )


## Examples

```r
data(mini_giotto_single_cell)

spatPlot2D(mini_giotto_single_cell,cell_color = 'cell_types', point_size = 5)
```


