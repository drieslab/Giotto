from calendar import c
import os

import anndata as ad
import bento as bt
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from geopandas import GeoDataFrame
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


def create_AnnData(trainscripts: pd.DataFrame, cell_shape: GeoDataFrame, nucleus_shape) -> AnnData:
    cell_seg = create_seg_df(cell_shape, x='x', y='y', cell_id='cell')
    nucleus_seg = create_seg_df(nucleus_shape, x='x', y='y', cell_id='nucleus')
    return bt.io.prepare(molecules=trainscripts, cell_seg=cell_seg, other_seg={'nucleus': nucleus_seg})  # type: ignore
