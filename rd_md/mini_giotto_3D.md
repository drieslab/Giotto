# `mini_giotto_3D`

mini Giotto object for spatial single-cell 3D data


## Description

Mini Giotto object created from the STARmap data.


## Format

An object of class `"giotto"` ; see [`createGiottoObject`](#creategiottoobject) .


## Usage

```r
data(mini_giotto_3D)
```


## References

Wang et al. (2018) Science
 ( [PubMed](https://pubmed.ncbi.nlm.nih.gov/29930089/) )


## Examples

```r
data(mini_giotto_3D)

spatPlot3D(mini_giotto_3D, cell_color = 'cell_types', point_size = 5)
```


