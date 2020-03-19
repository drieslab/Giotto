
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- This line is from RStudio -->

# Giotto

<!-- badges: start -->

<!-- badges: end -->

The Giotto package provides tools to process, analyze and visualize
**single-cell spatial expression** data. The underlying framework is
generalizable to virtually all currently available spatial datasets. We
recently demonstrated the general applicability on 10 different datasets
created by 9 different state-of-the-art spatial technologies, including
*in situ* hybridization (seqFISH+, merFISH, osmFISH), sequencing
(Slide-seq, Visium, STARmap) and imaging-based multiplexing/proteomics
(CyCIF, MIBI, CODEX). These technologies differ in terms of resolution
(single cell vs multiple cells), spatial dimension (2D vs 3D), molecular
modality (protein vs RNA), and throughput (number of cells and genes).

## Requirements

  - R (\>= 3.5.1)
  - Python (\>= 3.0)
  - Windows, MacOS, Linux

 

## Installation

See FAQs for additional information.

#### R installation

You can install (\~1-5 mins) Giotto with:

``` r
library(devtools) # if not installed: install.packages('devtools') 
library(remotes)  # if not installed: install.packages('remotes')
# to install the latest version (development)
remotes::install_github("RubD/Giotto")

# to install a specific release version (example)
remotes::install_github("RubD/Giotto@v0.1.3")
```

#### Python tools (optional)

This is necessary to run all available analyses, including Leiden /
Louvain clustering and to build and use the interactive visualization
tool. An alternative, but less flexible, R version for Louvain
clustering is also available. It is advisable to install everything
within a specific conda environment and specify the python path at the
beginning with createGiottoInstructions() or in the R function itself
when required.

Required python modules: pandas / igraph / networkx / leidenalg

pip installation one-liner:

``` bash
pip3 install pandas python-igraph networkx python-louvain leidenalg
```

If pip install does not work, try installing within a [conda
environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands):

``` bash
conda install -c anaconda pandas
conda install -c conda-forge python-igraph
conda install -c anaconda networkx
conda install -c conda-forge python-louvain
conda install -c conda-forge leidenalg
```

#### HMRF

See [**HMRF**](http://www.nature.com/articles/nbt.4260) installation
[instructions](http://spatial.rc.fas.harvard.edu/install.html).

#### Giotto Viewer

``` bash
pip3 install --user jsbeautifier
pip3 install --user giotto-viewer --no-cache #add --no-deps if do not wish to upgrade dependency
pip3 install --user smfish-image-processing --no-cache #add --no-deps if do not wish to upgrade dependency
```

 

## Examples

  - see <https://github.com/RubD/spatial-datasets> to find raw and
    pre-processed input data and Giotto scripts (in progress).
  - typical run time range for the different datasets on a personal
    computer is around 10\~45 mins.  
  - click on the image and try them out yourself.

[![Cortex](./inst/images/cortex_image_summary.png)](./inst/examples/mouse_cortex_svz/seqfish_cortex_Giotto_v0.1.4.md)
[![STARmap](./inst/images/starmap_cortex_image_summary.png)](./inst/examples/mouse_starmap_cortex/starmap_cortex_Giotto_v0.1.2_update.md)

[![osmFISH](./inst/images/osmFISH_SS_cortex_image_summary.png)](./inst/examples/mouse_osmFISH_SS_cortex/osmFISH_SS_cortex_Giotto_v0.1.2.md)
[![Spatial
Transcriptomics](./inst/images/SpatTx_OB_image_summary.png)](./inst/examples/mouse_SpatTx_OB/mouse_SpatTx_OB_Giotto_v0.1.2.md)

[![visium-kidney](./inst/images/visium_kidney_image_summary.png)](./inst/examples/mouse_visium_kidney/mouse_visium_kidney_Giotto_v0.1.3.md)
[![visium-brain](./inst/images/visium_brain_image_summary.png)](./inst/examples/mouse_visium_brain/mouse_visium_brain_Giotto_v0.1.3.md)

[![slideseq](./inst/images/slideseq_cerebellum_image_summary.png)](./inst/examples/mouse_slideseq_cerebellum/slideseq_cerebellum_Giotto_v0.1.2.md)
[![merFISH](./inst/images/merFISH_hypothalam_image_summary.png)](./inst/examples/mouse_hypothalamic_preoptic_region/merfish_hypo_Giotto_v0.1.2_update.md)

\[![MIBI-TOF](./inst/images/MIBI_tumor_image_summary.png)\]

 

## References

  - Dries, R. et al. Giotto, a pipeline for integrative analysis and
    visualization of single-cell spatial transcriptomic data. bioRxiv
    701680 (2019). <doi:10.1101/701680>
    [link](https://www.biorxiv.org/content/10.1101/701680v1)

  - Eng, C.-H. L. et al. Transcriptome-scale super-resolved imaging in
    tissues by RNA seqFISH+. Nature 1 (2019).
    <doi:10.1038/s41586-019-1049-y>
    [link](https://www.nature.com/articles/s41586-019-1049-y)

  - Zhu, Q., Shah, S., Dries, R., Cai, L. & Yuan, G.-C. Identification
    of spatially associated subpopulations by combining scRNAseq and
    sequential fluorescence in situ hybridization data. Nature
    Biotechnology (2018). <doi:10.1038/nbt.4260>
    [link](https://www.nature.com/articles/nbt.4260)
