# `getSpatialDataset`

getSpatialDataset


## Description

This package will automatically download the spatial locations and
 expression matrix for the chosen dataset. These files are already in the right format
 to create a Giotto object. If wget is installed on your machine, you can add
 'method = wget' to the parameters to download files faster.


## Usage

```r
getSpatialDataset(
  dataset = c("ST_OB1", "ST_OB2", "codex_spleen", "cycif_PDAC", "starmap_3D_cortex",
    "osmfish_SS_cortex", "merfish_preoptic", "seqfish_SS_cortex", "seqfish_OB",
    "slideseq_cerebellum", "Human_PCa_scRNAseq", "Mouse_brain_scRNAseq"),
  directory = getwd(),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`dataset`     |     dataset to download
`directory`     |     directory to save the data to
`list()`     |     additional parameters to [`download.file`](#download.file)


