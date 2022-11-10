import anndata as ad
import pandas as pd
import numpy as np
import scipy

### Imports and safeguards
def read_anndata_from_path(ad_path = None):
    if ad_path is None:
        print("Please provide a path to a .h5ad file.")
        assert(False)
    try:
        adata = ad.read(ad_path)
    except (FileNotFoundError):
        print(f"File {ad_path} was not found.")
        print("Please try again.")
        raise(FileNotFoundError)
    finally:
        adata.file.close()
    return adata

def ad_guard(adata):
    if type(adata) is not type(ad.AnnData()): 
        print("Please provide a valid AnnData object.")
        raise(TypeError)

### Expression info
def extract_expression(adata = None):
    ad_guard(adata)
    # anndata stores expression as cells by feats
    # giotto stores expression as feats by cells
    expr = adata.X.transpose()
    return expr

### IDs
def extract_cell_IDs(adata = None):
    ad_guard(adata)
    cell_IDs = [i for i in adata.obs_names]
    return cell_IDs

def extract_feat_IDs(adata = None):
    ad_guard(adata)
    feat_IDs = [i for i in adata.var_names]
    return feat_IDs

### Dimension reductions
def extract_pca(adata = None):
    """
    Extracts PCA information from AnnData object
        - PCA array
        - Loadings
        - Eigenvalues of covariance matrix
    
    INPUT: AnnData object

    OUTPUT: Dictionary containing PCA information
    """


    ad_guard(adata)
    o_keys = adata.obsm_keys()
    v_keys = adata.varm_keys()
    u_keys = adata.uns['pca'].keys()

    pca = dict()

    for ok in o_keys:
        if "pca" in ok or "PCA" in ok:
            pca['pca'] = adata.obsm[ok]
    for vk in v_keys:
        if "PC" in vk:
            pca['loadings'] = adata.varm[vk]
    for uk in u_keys:
        if "variance" == uk:
            pca['eigenvalues'] = adata.uns['pca'][uk]
    
    return pca
    
def extract_umap(adata = None):
    ad_guard(adata)
    o_keys = adata.obsm_keys()

    umap = None

    for ok in o_keys:
        if "umap" in ok or "UMAP" in ok:
            umap = adata.obsm[ok]

    return umap  

def extract_tsne(adata = None):
    ad_guard(adata)
    o_keys = adata.obsm_keys()

    tsne = None

    for ok in o_keys:
        if "tsne" in ok or "tsne" in ok:
            tsne = adata.obsm[ok]

    return tsne  

### Spatial
def parse_obsm_for_spat_locs(adata = None):
    """
    Parses anndata.AnnData.obsm for spatial location information
    
    NOTE: obsm keyword "spatial" must be present for this function to work. The word
    'spatial' MUST at least appear within the keyword 
    
    i.e. 'spatial_location' works, 'spat_location' does not work
        
    INPUT: AnnData object

    OUTPUT: numpy array containing centroid information (x,y(,z))
    """
    ad_guard(adata)
    cID = np.array(extract_cell_IDs(adata))
    spat_locs = None
    try:
        spat_locs = adata.obsm["spatial"]
    except (KeyError):
        spat_key = [i for i in adata.obsm if 'spatial' in i][0]
        spat_locs = adata.obsm[spat_key]
    
    if spat_locs is None:
        print("Spatial locations were not found.")
        raise(KeyError)
    else:
        print("Spatial locations found.")
    
    cID = np.array(cID).reshape(len(cID),1)
    spat_locs = np.concatenate((spat_locs,cID), axis = 1)
    num_col = spat_locs.shape[1]
    
    colnames = ["sdimx","sdimy","sdimz","cell_ID"]
    conv = {"sdimx":float, "sdimy":float,"sdimz":float,"cell_ID":str}
    
    if num_col > 3:
        spat_locs = pd.DataFrame(spat_locs, columns = colnames)
    else:
        del colnames[2]
        del conv['sdimz']
        spat_locs = pd.DataFrame(spat_locs, columns = colnames)
    
    spat_locs = spat_locs.astype(conv)
    return spat_locs

### Metadata
def extract_cell_metadata(adata = None):
    ad_guard(adata)
    cell_metadata = adata.obs.reset_index()
    cell_metadata = cell_metadata.rename(columns={"index":"cell_ID"})
    return cell_metadata

def extract_feat_metadata(adata = None):
    ad_guard(adata)
    feat_metadata = adata.var.reset_index()
    feat_metadata = feat_metadata.rename(columns={"index":"feat_ID"})
    return feat_metadata

### Alternative expression data
def extract_layer_names(adata = None):
    ad_guard(adata)
    layer_names = [i for i in adata.layers]
    return layer_names

def extract_layered_data(adata = None, layer_name = None):
    ad_guard(adata)
    layer_names = [i for i in adata.layers]
    if layer_name not in layer_names:
        print(f"Invalid Key, {layer_name}, for adata.layers")
        raise(KeyError)
    target_layer = adata.layers[layer_name]
    if type(target_layer) == scipy.sparse.csr.csr_matrix:
        target_layer = target_layer.T
    else:
        target_layer = pd.DataFrame(target_layer)
    return target_layer

def unstruct(adata = None):
    """
    Extracts unstructured data from AnnData object

    INPUT: AnnData object

    OUTPUT: Dictionary of dictionaries containing unstructured metadata
    """
    ad_guard(adata)
    unstructured = adata.uns

    return unstructured
