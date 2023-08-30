from typing import List, Optional

import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csc_matrix
from shapely.geometry import MultiPolygon, Point, Polygon


def create_seg_df(vertices_df: pd.DataFrame, x: str = 'x', y: str = 'y', cell_id: str = 'cell_id') -> pd.DataFrame:
    """
    Create a dataframe with cell_id and geometry columns from a dataframe with vertices
    :param vertices_df: a dataframe with columns: cell_id, x, y
    :param x: the column name of x coordinates
    :param y: the column name of y coordinates
    :param cell_id: the column name of cell id
    :param bounds: the bounds of the area (minx, miny, maxx, maxy)
    :return: a dataframe with cell_id and geometry columns
    """
    # --- create polygons ---
    polygons = vertices_df.groupby(cell_id).apply(lambda group: Polygon(zip(group[x], group[y])))  # type: ignore
    seg_df = gpd.GeoDataFrame(polygons, columns=['geometry'])
    # --- correct invalid polygons ---
    corrected_seg_df = seg_df.copy(deep=True)
    for i in range(seg_df.shape[0]):
        if not seg_df.iloc[i, 0].is_valid:  # type: ignore
            corrected_seg_df.iloc[i, 0] = seg_df.iloc[i, 0].buffer(0)  # type: ignore
        if isinstance(corrected_seg_df.iloc[i, 0], MultiPolygon):
            corrected_seg_df.iloc[i, 0] = max(corrected_seg_df.iloc[i, 0], key=lambda x: x.area)  # type: ignore
    return corrected_seg_df


def valid_transcript(row, adata: AnnData, legal_cells: List[str]) -> bool:
    """
    Check if the transcript is valid
    :param adata: anndata object
    :param row: a row of adata.uns['points']
    :return: True if the transcript is valid
    """
    return str(row['cell']) in legal_cells and adata.obs.loc[str(row['cell']), 'cell_shape'].contains(  # type: ignore
        Point(row['x'], row['y']))


def adata_X_handle(adata: AnnData) -> None:
    """
    Convert adata.X from sparse matrix to dense matrix
    :param adata: anndata object
    :return: None
    """
    if isinstance(adata.X, csc_matrix):  # in case of duplicate handling
        adata.X = np.array(adata.X.todense())


def adata_var_handle(adata: AnnData) -> None:
    """
    Set gene as index of adata.var
    :param adata: anndata object
    :return: None
    """
    if adata.var.shape[1] == 1:  # in case of duplicate handling
        adata.var.set_index('gene', inplace=True)


def adata_obs_handle(adata: AnnData,
                     shape_names: Optional[List[str]] = None,
                     bounds: Optional[List[float]] = None) -> None:
    """
    Handle adata.obs to make sure the data type is correct
    :param adata: anndata object
    :param shape_names: the names of shapes that need to handle
    :param bounds: the bounds of the area (minx, miny, maxx, maxy)
    :return: None
    """
    # --- shapes need to handle ---
    if shape_names is None:
        shape_names = ['cell_shape', 'nucleus_shape']  # we may have other shapes (mitochondria etc) in the future
    # --- create GEODataFrame ---
    for shape in shape_names:
        if shape not in adata.obs.columns:
            seg_df = create_seg_df(vertices_df=adata.uns[shape])
            adata.obs[shape] = seg_df['geometry']
            del adata.uns[shape]
    # --- remove cells outside the bounds ---
    if bounds is not None:
        area_shape = Polygon([(bounds[0], bounds[1]), (bounds[0], bounds[3]), (bounds[2], bounds[3]),
                              (bounds[2], bounds[1])])
        adata.uns['legal_cells'] = adata.obs['cell_shape'].apply(lambda x: area_shape.contains(x))
    else:
        adata.uns['legal_cells'] = np.ones(adata.obs.shape[0], dtype=bool)
    # --- add batch column ---
    if 'batch' not in adata.obs.columns:
        adata.obs['batch'] = 0


def adata_uns_points_handle(adata: AnnData) -> None:
    """
    Handle adata.uns['points'] to make sure the data type is correct
    :param adata: anndata object
    :return: None
    """
    # --- correct the value of nucleus ---
    if adata.uns['points']['nucleus'].min() == 0:  # in case of duplicate handling
        points = adata.uns['points']
        points['nucleus'] = points['cell'].values * points['nucleus'].values
        points.loc[points['nucleus'] == 0, 'nucleus'] = -1
    # --- remove unused transcripts ---
    legal_cells = adata.obs.index[adata.uns['legal_cells']]
    legal_transcripts = adata.uns['points'].apply(valid_transcript, axis=1, args=(adata, legal_cells))
    adata.uns['points'] = adata.uns['points'][legal_transcripts]
    # --- setup data type ---
    dtypes = {'x': 'float', 'y': 'float', 'gene': 'str', 'cell': 'str', 'nucleus': 'str', 'batch': 'str'}
    adata.uns['points'] = adata.uns['points'].astype(dtypes)
    adata.uns['points']['gene'] = adata.uns['points']['gene'].astype('category')
    adata.uns['points']['cell'] = adata.uns['points']['cell'].astype('category')
    adata.uns['points']['nucleus'] = adata.uns['points']['nucleus'].astype('category')
    adata.uns['points']['batch'] = adata.uns['points']['batch'].astype('category')
    # --- add legal_genes in adata.uns ---
    adata.uns['legal_genes'] = [x in adata.uns['points']['gene'].tolist() for x in adata.var.index]


def adata_handle(adata: AnnData, bounds: Optional[List[float]] = None) -> AnnData:
    """
    Handle anndata object to make sure the data type is correct
    Process:
        1. processing adata.X: format
        2. processing adata.var: format
        3. processing adata.obs
            - create GEODataFrame & delete adata.uns['cell_shape'] and adata.uns['nucleus_shape']
            - remove cells outside the bounds & add legal_cells in adata.uns
            - add batch column
        4. processing adata.uns['points']
            - correct the value of nucleus
            - remove unused transcripts
            - setup data type
            - add legal_genes in adata.uns
        5. subset adata
    :param adata: AnnData object
    :return: Filtered Anndata object
    """
    adata_X_handle(adata)
    adata_var_handle(adata)
    adata_obs_handle(adata, bounds=bounds)
    adata_uns_points_handle(adata)

    return adata[adata.uns['legal_cells'], adata.uns['legal_genes']]


