import bento as bt
import geopandas as gpd
import pandas as pd
from anndata import AnnData
from shapely.geometry import MultiPolygon, Polygon

from log import debug, warning, info

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


def add_batch(adata: AnnData, cell_shape: pd.DataFrame):
    """
    Add batch information to an AnnData object
        If cell_seg have batch information, add batch information to adata, else all batch will be set to 0
    :param adata: an AnnData object
    :param cell_seg: the name of the cell segmentation
    :return: an AnnData object with batch information
    """
    if 'batch' in cell_shape.columns:
        info('Batch information found in cell_shape, adding batch information to adata')
        adata.obs['batch'] = [cell_shape.loc[cell_shape['cell_id']==cell,'batch'].values[0] for cell in adata.obs_names] # type: ignore
        adata.uns['points']['batch'] = [cell_shape.loc[cell_shape['cell_id']==cell,'batch'].values[0] for cell in adata.uns['points']['cell']] # type: ignore
    else:  # Interim measures, batch information may not transfered to cell_shape
        warning('Batch information not found in cell_shape, all batch will be set to 0')
        adata.obs['batch'] = 0
        adata.uns['points']['batch'] = 0
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    adata.uns['points']['batch'] = adata.uns['points']['batch'].astype('category')


def create_AnnData(trainscripts, cell_shape, nucleus_shape) -> AnnData:
    # --- processing input ---
    trainscripts = pd.DataFrame(trainscripts)
    cell_shape = pd.DataFrame(cell_shape)
    cell_shape['cell_id'] = cell_shape['cell_id'].astype('category')
    if 'batch' in cell_shape.columns:
        cell_shape['batch'] = cell_shape['batch'].astype('category')
    nucleus_shape = pd.DataFrame(nucleus_shape)
    nucleus_shape['cell_id'] = nucleus_shape['cell_id'].astype('category')

    # --- create shape ---
    cell_seg = create_seg_df(cell_shape, x='x', y='y', cell_id='cell_id')
    nucleus_seg = create_seg_df(nucleus_shape, x='x', y='y', cell_id='cell_id')
    if cell_seg.shape[0] > 500:
        warning('cell_seg has more than 500 cells, processing may take a long time.')

    # --- filter cells ---
    # Let Giotto perform the filtering
    legal_cells = pd.Series([True] * cell_seg.shape[0])

    # --- create AnnData ---
    adata: AnnData = bt.io.prepare(molecules=trainscripts, cell_seg=cell_seg, other_seg={'nucleus': nucleus_seg})  # type: ignore
    add_batch(adata, cell_shape)

    # --- filter genes ---
    # Interim measures
    # subsetGiottoLocs don't give wanted result, so we use a workaround here
    legal_genes = adata.var_names.isin(set(adata.uns['points']['gene'].values))
    filtered_adata = adata[legal_cells,legal_genes] # type: ignore
    filtered_adata.uns['points']['gene'] = filtered_adata.uns['points']['gene'].cat.remove_unused_categories()

    # Interim measures
    # subsetted adata (_is_view == True) don't have _X property, which will cause unexpect error for adata.__sizeof__
    # when create object in R, reticulate will call sys.getsizeof to get the size of the object, which will call adata.__sizeof__
    # https://github.com/rstudio/reticulate/issues/1332
    # https://github.com/rstudio/rstudio/issues/13491
    # repoted to AnnData here: https://github.com/scverse/anndata/issues/1127
    filtered_adata._X = filtered_adata.X 
    return filtered_adata

    # return adata
