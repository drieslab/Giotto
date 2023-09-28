from typing import List, Optional, Union, Iterable
import bento as bt
from anndata import AnnData
import emoji
from bento._utils import track
from bento.tools._colocation import _colocation_tensor
from bento.tools import decompose
from kneed import KneeLocator
from sklearn.utils import resample
from minisom import MiniSom
from tqdm.auto import tqdm
import numpy as np
from scipy.sparse import csr_matrix
import rasterio
import shapely
import geopandas as gpd
import pandas as pd
from bento.geometry import sindex_points

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
            fig, ax = plt.subplots(figsize=(6, 4))
            sns.lineplot(data=errors, x="rank", y="rmse", ci=95, marker="o", ax=ax)  # type: ignore
            plt.axvline(kl.knee, linestyle="--")
            plt.savefig(fname)
            plt.close(fig)
            info(f"Saved to {fname}")

    adata.uns["factors"] = factors
    adata.uns["factors_error"] = errors

    print(emoji.emojize(":heavy_check_mark: Done."))
    return adata if copy else None


@track
def fluxmap(
    data: AnnData,
    fname: str,
    n_clusters: Union[Iterable[int], int] = range(2, 9),
    num_iterations: int = 1000,
    train_size: float = 0.2,
    res: float = 0.1,
    random_state: int = 11,
    plot_error: bool = True,
    copy: bool = False,
):
    """Cluster flux embeddings using self-organizing maps (SOMs) and vectorize clusters as Polygon shapes.

    Parameters
    ----------
    data : AnnData
        Spatial formatted AnnData object.
    n_clusters : int or list
        Number of clusters to use. If list, will pick best number of clusters
        using the elbow heuristic evaluated on the quantization error.
    num_iterations : int
        Number of iterations to use for SOM training.
    train_size : float
        Fraction of cells to use for SOM training. Default 0.2.
    res : float
        Resolution used for rendering embedding. Default 0.05.
    random_state : int
        Random state to use for SOM training. Default 11.
    plot_error : bool
        Whether to plot quantization error. Default True.
    copy : bool
        Whether to return a copy the AnnData object. Default False.

    Returns
    -------
    adata : AnnData
        .uns["cell_raster"] : DataFrame
            Adds "fluxmap" column denoting cluster membership.
        .uns["points"] : DataFrame
            Adds "fluxmap#" columns for each cluster.
        .obs : GeoSeries
            Adds "fluxmap#_shape" columns for each cluster rendered as (Multi)Polygon shapes.
    """
    adata = data.copy() if copy else data

    # Check if flux embedding has been computed
    if "flux_embed" not in adata.uns:
        raise ValueError("Flux embedding has not been computed. Run `bento.tl.flux()` first.")

    flux_embed = adata.uns["flux_embed"]
    raster_points = adata.uns["cell_raster"]

    if isinstance(n_clusters, int):
        n_clusters = [n_clusters]

    if isinstance(n_clusters, range):
        n_clusters = list(n_clusters)

    # Subsample flux embeddings for faster training
    if train_size > 1:
        raise ValueError("train_size must be less than 1.")
    if train_size == 1:
        flux_train = flux_embed
    if train_size < 1:
        flux_train = resample(
            flux_embed,
            n_samples=int(train_size * flux_embed.shape[0]),
            random_state=random_state,
        )

    # Perform SOM clustering over n_clusters range and pick best number of clusters using elbow heuristic
    pbar = tqdm(total=4)
    pbar.set_description(emoji.emojize(f"Optimizing # of clusters"))
    som_models = {}
    quantization_errors = []
    for k in tqdm(n_clusters, leave=False):
        som = MiniSom(1, k, flux_train.shape[1], random_seed=random_state)  # type: ignore
        som.random_weights_init(flux_train)  # type: ignore
        som.train(flux_train, num_iterations, random_order=False, verbose=False)  # type: ignore
        som_models[k] = som
        quantization_errors.append(som.quantization_error(flux_embed))

    # Use kneed to find elbow
    if len(n_clusters) > 1:  # type: ignore
        kl = KneeLocator(n_clusters, quantization_errors, curve="convex", direction="decreasing")
        if kl.elbow is None:
            warning(
                'No elbow found, please extend the n_clusters range.\nCurrent n_clusters range: [{n_clusters[0]},{n_clusters[-1]}]'
            )
            return adata if copy else None
        else:
            info(f'Elbow found at {kl.elbow}')
            best_k = kl.elbow
            fig, ax = plt.subplots(figsize=(6, 4))
            sns.lineplot(x=n_clusters, y=quantization_errors, ci=95, marker="o", ax=ax)  # type: ignore
            plt.axvline(kl.elbow, linestyle="--")
            plt.savefig(fname)
            plt.close(fig)
            info(f"Saved to {fname}")

    else:
        best_k = n_clusters[0]  # type: ignore
    pbar.update()

    # Use best k to assign each sample to a cluster
    pbar.set_description(f"Assigning to {best_k} clusters")
    som = som_models[best_k]
    winner_coordinates = np.array([som.winner(x) for x in flux_embed]).T

    # Indices start at 0, so add 1
    qnt_index = np.ravel_multi_index(winner_coordinates, (1, best_k)) + 1  # type: ignore
    raster_points["fluxmap"] = qnt_index
    adata.uns["cell_raster"] = raster_points.copy()

    pbar.update()

    # Vectorize polygons in each cell
    pbar.set_description(emoji.emojize("Vectorizing domains"))
    cells = raster_points["cell"].unique().tolist()
    # Scale down to render resolution
    # raster_points[["x", "y"]] = raster_points[["x", "y"]] * res

    # Cast to int
    raster_points[["x", "y", "fluxmap"]] = raster_points[["x", "y", "fluxmap"]].astype(int)

    rpoints_grouped = raster_points.groupby("cell")
    fluxmap_df = dict()
    for cell in tqdm(cells, leave=False):
        rpoints = rpoints_grouped.get_group(cell)

        # Fill in image at each point xy with fluxmap value by casting to dense matrix
        image = (csr_matrix((
            rpoints["fluxmap"],
            (
                (rpoints["y"] * res).astype(int),
                (rpoints["x"] * res).astype(int),
            ),
        )).todense().astype("int16"))

        # Find all the contours
        contours = rasterio.features.shapes(image)  # type: ignore
        polygons = np.array([(shapely.geometry.shape(p), v) for p, v in contours])  # type: ignore
        shapes = gpd.GeoDataFrame(
            polygons[:, 1],
            geometry=gpd.GeoSeries(polygons[:, 0]).T,
            columns=["fluxmap"],
        )  # type: ignore

        # Remove background shape
        shapes["fluxmap"] = shapes["fluxmap"].astype(int)  # type: ignore
        shapes = shapes[shapes["fluxmap"] != 0]

        # Group same fields as MultiPolygons
        shapes = shapes.dissolve("fluxmap")["geometry"]  # type: ignore

        fluxmap_df[cell] = shapes

    fluxmap_df = pd.DataFrame.from_dict(fluxmap_df).T
    fluxmap_df.columns = "fluxmap" + fluxmap_df.columns.astype(str) + "_shape"

    # Upscale to match original resolution
    fluxmap_df = fluxmap_df.apply(lambda col: gpd.GeoSeries(col).scale(xfact=1 / res, yfact=1 / res, origin=(0, 0)))
    pbar.update()

    pbar.set_description("Saving")
    old_cols = adata.obs.columns[adata.obs.columns.str.startswith("fluxmap")]
    adata.obs = adata.obs.drop(old_cols, axis=1, errors="ignore")

    adata.obs[fluxmap_df.columns] = fluxmap_df.reindex(adata.obs_names)

    old_cols = adata.uns["points"].columns[adata.uns["points"].columns.str.startswith("fluxmap")]
    adata.uns["points"] = adata.uns["points"].drop(old_cols, axis=1)

    # TODO SLOW
    sindex_points(adata, "points", fluxmap_df.columns.tolist())
    pbar.update()
    pbar.set_description("Done")
    pbar.close()

    return adata if copy else None


# ---------------------------------
# bento wrapper functions
# ---------------------------------


def analysis_shape_features(adata: AnnData, feature_names: Optional[List[str]] = None) -> None:
    """
    Examples
    --------
    >>> analysis_shape_features(adata=bento_adata)
    """
    if feature_names is None:
        feature_names = list(bt.tl.list_shape_features().keys())
    bt.tl.obs_stats(adata, feature_names=feature_names)


def plot_shape_features_analysis_results(adata: AnnData, fname: str):
    """
    Examples
    --------
    >>> plot_shape_features_analysis_results(adata=bento_adata, fname='test_shape_features.pdf')
    """
    bt.pl.shapes(adata, fname=fname)


def analysis_points_features(adata: AnnData,
                             shapes_names: Optional[List[str]] = None,
                             feature_names: Optional[List[str]] = None) -> None:
    """
    Examples
    --------
    >>> analysis_points_features(adata=bento_adata)
    """
    if shapes_names is None:
        shapes_names = ["cell_shape", "nucleus_shape"]
    if feature_names is None:
        feature_names = list(bt.tl.list_point_features().keys())
    bt.tl.analyze_points(adata, shape_names=shapes_names, feature_names=feature_names, groupby='gene')


def plot_points_features_analysis_results(adata: AnnData, fname: str) -> None:
    """
    Examples
    --------
    >>> plot_points_features_analysis_results(adata=bento_adata, fname='test_points_features.pdf')
    """
    bt.pl.points(adata, fname=fname)


def analysis_RNAflux(adata: AnnData) -> None:
    """
    Examples
    --------
    >>> analysis_RNAflux(adata=bento_adata)
    """
    bt.tl.flux(adata)


def plot_RNAflux_analysis_results(adata: AnnData, fname: str) -> None:
    """
    Examples
    --------
    >>> plot_RNAflux_analysis_results(adata=bento_adata, fname='test_RNAflux.pdf')
    """
    bt.pl.flux(adata, fname=fname)


def analysis_RNAfluxmap(adata: AnnData, fname: str, n_clusters: Union[Iterable[int], int] = range(2, 9)) -> None:
    """
    Examples
    --------
    >>> analysis_RNAfluxmap(adata=bento_adata, fname='test_RNAfluxmap_elbow_pos.png', n_clusters=seq(20))
    """
    if fname is None:
        fname = 'fluxmap_elbow.pdf'
    fluxmap(adata, fname, n_clusters=n_clusters)


def plot_RNAfluxmap_analysis_results(adata: AnnData, fname: str) -> None:
    """
    Examples
    --------
    >>> plot_RNAfluxmap_analysis_results(adata=bento_adata, fname='test_RNAfluxmap.pdf')
    """
    bt.pl.fluxmap(adata, fname=fname)


def analysis_rna_forest(adata: AnnData) -> None:
    """
    Examples
    --------
    >>> analysis_rna_forest(adata=bento_adata)
    """
    bt.tl.lp(adata)
    bt.tl.lp_stats(adata)


def plot_rna_forest_analysis_results(adata: AnnData, fname1: str, fname2: str) -> None:
    """
    Examples
    --------
    >>> plot_rna_forest_analysis_results(adata=bento_adata, fname1='test_rna_forest_radvis.pdf', fname2='test_rna_forest_upset.pdf')
    """
    bt.pl.lp_genes(adata, fname=fname1)
    bt.pl.lp_dist(adata, fname=fname2)


def analysis_colocalization(adata: AnnData, fname: str, ranks: Optional[List[int]] = None) -> None:
    """
    Examples
    --------
    >>> analysis_colocalization(adata=bento_adata, fname='test_colocalization_knee_pos.pdf')
    """
    if ranks is None:
        ranks = list(range(1, 6))

    # Cytoplasm = cell - nucleus
    adata.obs["cytoplasm_shape"] = bt.geo.get_shape(adata, "cell_shape") - bt.geo.get_shape(adata, "nucleus_shape")

    # Create point index
    adata.uns["points"]["cytoplasm"] = (adata.uns["points"]["nucleus"].astype(int) < 0).astype(int)

    bt.tl.coloc_quotient(adata, shapes=["cytoplasm_shape", "nucleus_shape"])

    colocation(adata, ranks=ranks, fname=fname)


def plot_colocalization_analysis_results(adata: AnnData, fname: str, rank: int) -> None:
    """
    Examples
    --------
    >>> plot_colocalization_analysis_results(adata=bento_adata, fname='test_colocalization.pdf', rank=2)
    """
    bt.pl.colocation(adata, rank=rank, fname=fname)


def python_session_info():
    import session_info
    session_info.show(html=False, dependencies=True)
