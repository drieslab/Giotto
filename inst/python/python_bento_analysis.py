from typing import List, Optional
import bento as bt
from anndata import AnnData
import emoji
from bento._utils import track
from bento.tools._colocation import _colocation_tensor
from bento.tools import decompose
from kneed import KneeLocator

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from log import warning, info


# ---------------------------------
# modified bento and dependencies functions/classes
# ---------------------------------
@track
def colocation(
    data: AnnData,
    ranks: List[int],
    fname: str,
    iterations: int = 3,
    plot_error: bool = True,
    copy: bool = False,
):
    """Decompose a tensor of pairwise colocalization quotients into signatures.

    Parameters
    ----------
    adata : AnnData
        Spatial formatted AnnData object.
    ranks : list
        List of ranks to decompose the tensor.
    iterations : int
        Number of iterations to run the decomposition.
    plot_error : bool
        Whether to plot the error of the decomposition.
    copy : bool
        Whether to return a copy of the AnnData object. Default False.
    Returns
    -------
    adata : AnnData
        .uns['factors']: Decomposed tensor factors.
        .uns['factors_error']: Decomposition error.
    """
    adata = data.copy() if copy else data

    print("Preparing tensor...")
    _colocation_tensor(adata, copy=copy)

    tensor = adata.uns["tensor"]

    print(emoji.emojize(":running: Decomposing tensor..."))
    factors, errors = decompose(tensor, ranks, iterations=iterations)

    if plot_error and errors.shape[0] > 1:
        kl = KneeLocator(errors["rank"], errors["rmse"], direction="decreasing", curve="convex")
        if kl.knee is None:
            warning('No knee found, please extend the ranks range.\nCurrent ranks range: [{ranks[0]},{ranks[-1]}]')
        else:
            info(f'Knee found at rank {kl.knee}')
            sns.lineplot(data=errors, x="rank", y="rmse", ci=95, marker="o")  # type: ignore
            plt.axvline(kl.knee, linestyle="--")
            plt.savefig(fname)
            info(f"Saved to {fname}")

    adata.uns["factors"] = factors
    adata.uns["factors_error"] = errors

    print(emoji.emojize(":heavy_check_mark: Done."))
    return adata if copy else None


# ---------------------------------
# bento wrapper functions
# ---------------------------------


def analysis_shape_features(adata: AnnData, feature_names: Optional[List[str]] = None) -> None:
    if feature_names is None:
        feature_names = list(bt.tl.list_shape_features().keys())
    bt.tl.obs_stats(adata, feature_names=feature_names)


def plot_shape_features_analysis_results(adata: AnnData, fname: str):
    bt.pl.shapes(adata, fname=fname)


def analysis_points_features(adata: AnnData,
                             shapes_names: Optional[List[str]] = None,
                             feature_names: Optional[List[str]] = None) -> None:
    if shapes_names is None:
        shapes_names = ["cell_shape", "nucleus_shape"]
    if feature_names is None:
        feature_names = list(bt.tl.list_point_features().keys())
    bt.tl.analyze_points(adata, shape_names=shapes_names, feature_names=feature_names, groupby='gene')


def plot_points_features_analysis_results(adata: AnnData, fname: str) -> None:
    bt.pl.points(adata, fname=fname)


def analysis_rna_forest(adata: AnnData) -> None:
    bt.tl.lp(adata)
    bt.tl.lp_stats(adata)


def plot_rna_forest_analysis_results(adata: AnnData, fname1: str, fname2: str) -> None:
    bt.pl.lp_genes(adata, fname=fname1)
    bt.pl.lp_dist(adata, fname=fname2)


def analysis_colocalization(adata: AnnData, fname: str, ranks: Optional[List[int]] = None) -> None:
    if ranks is None:
        ranks = list(range(1, 6))

    # Cytoplasm = cell - nucleus
    adata.obs["cytoplasm_shape"] = bt.geo.get_shape(adata, "cell_shape") - bt.geo.get_shape(adata, "nucleus_shape")

    # Create point index
    adata.uns["points"]["cytoplasm"] = (adata.uns["points"]["nucleus"].astype(int) < 0).astype(int)

    bt.tl.coloc_quotient(adata, shapes=["cytoplasm_shape", "nucleus_shape"])

    colocation(adata, ranks=ranks, fname=fname)


def plot_colocalization_analysis_results(adata: AnnData, fname: str, rank: int) -> None:
    bt.pl.colocation(adata, rank=rank, fname=fname)


def chekc_genes_number(adata: AnnData) -> None:
    print(f'adata shape: {adata.shape}')
    print(f'adata points genes: {len(adata.uns["points"]["gene"].unique())}')
    print(f'adata cell_gene_features genes: {len(adata.uns["cell_gene_features"]["gene"].unique())}')
    diff_set = set(adata.uns["points"]["gene"]) - set(adata.uns["cell_gene_features"]["gene"])
    print(diff_set)
    for g in diff_set:
        print(f'{g}')
        print(adata.uns['points'][adata.uns['points']['gene'] == g])
