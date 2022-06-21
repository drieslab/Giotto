# `installGiottoEnvironment`

installGiottoEnvironment


## Description

Installs a giotto environment


## Usage

```r
installGiottoEnvironment(
  packages_to_install = c("pandas==1.1.5", "networkx==2.6.3", "python-igraph==0.9.6",
    "leidenalg==0.8.7", "python-louvain==0.15", "python.app==2", "scikit-learn==0.24.2"),
  force_miniconda = FALSE,
  force_environment = FALSE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`packages_to_install`     |     all python modules (packages) that should be installed for Giotto to work
`force_miniconda`     |     force reinstallation of miniconda
`force_environment`     |     force reinstallation of the giotto environment
`verbose`     |     be verbose


## Details

This function will install a local giotto environment using
 the miniconda system as implemented by reticulate. Once this giotto environment is
 installed it will be automatically detected when you run the Giotto toolbox. If you want to use
 your own python path then you can set the python_path in the [`createGiottoInstructions`](#creategiottoinstructions) 
 and provide the instructions to the [`createGiottoObject`](#creategiottoobject) function.


## Value

installs a giotto environment using the reticulate miniconda system


## Examples

```r
# this command will install r-miniconda
# and a giotto environment with all necessary python modules
installGiottoEnvironment()

The following python modules will be installed:
- pandas==1.1.5
- networkx==2.6.3
- python-igraph==0.9.6
- leidenalg==0.8.7
- python-louvain==0.15
- python.app==2 [mac OSX only]
- scikit-learn==0.24.2
```


