import anndata as ad
import pandas as pd
import numpy as np
import scipy
import warnings

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
    u_keys = None

    pca = dict()

    for ok in o_keys:
        if "pca" in ok or "PCA" in ok:
            pca['pca'] = adata.obsm[ok]
            u_keys = adata.uns['pca'].keys()
    for vk in v_keys:
        if "PC" in vk:
            pca['loadings'] = adata.varm[vk]
    if type(u_keys) is not type(None):
        for uk in u_keys:
            if "variance" == uk:
                pca['eigenvalues'] = adata.uns['pca'][uk]
    
    if(len(pca)) == 0:
        pca = None

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
    spat_key = None
    
    try:
        spat_locs = adata.obsm["spatial"]
    except (KeyError):
        spat_keys = [i for i in adata.obsm if 'spatial' in i]
        if len(spat_keys) > 0:
            spat_key = spat_keys[0]
            spat_locs = adata.obsm[spat_key]

    if spat_locs is None:
        err_mess = '''Spatial locations were not found. If spatial locations should have been found,
        please modify the anndata object to include a keyword-value pair within the obsm slot,
        in which the keyword contains the phrase "spatial" and the value corresponds to the spatial locations.\n
        In the Giotto Object resulting from this conversion, dummy locations will be used.'''
        print(err_mess)
        spat_locs = None
        return spat_locs
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
    layer_names = None
    if len(adata.layers) > 0:
        layer_names = [i for i in adata.layers]
    return layer_names

def extract_layered_data(adata = None, layer_name = None):
    ad_guard(adata)
    layer_names = [i for i in adata.layers]
    if layer_name not in layer_names:
        print(f"Invalid Key, {layer_name}, for adata.layers")
        raise(KeyError)
    target_layer = adata.layers[layer_name]
    if type(target_layer) == scipy.sparse.csr_matrix:
        target_layer = target_layer.T
    elif type(target_layer) == scipy.sparse.csr.csr_matrix:
        target_layer = target_layer.T
    else:
        target_layer = pd.DataFrame(target_layer)
    return target_layer

### Nearest Network

def find_NN_keys(adata = None, key_added = None):
    # No need to add to globals
    nn_key_list = []

    if key_added and key_added.casefold() != "spatial":
        map_keys = adata.uns[key_added].keys()
        for i in map_keys:
            #if type(adata.uns[key_added][i]) == type(dict()): continue
            nn_key_list.append(adata.uns[key_added][i])
    elif key_added and key_added.casefold() == "spatial":
        s1 = "String 'spatial' cannot be used as n_key_added to retrieve a Nearest Neighbor Network. "
        s2 = "This results from conflicting keys for nearest neighbor and spatial networks. "
        s3 = "\nSee defaults here:\nhttps://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html\nhttps://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html"
        msg = s1+ s2 + s3
        warnings.warn(msg)
    else: 
        param_keys = list(adata.uns.keys())
        for pk in param_keys:
            if "neighbors" in pk and "spatial" not in pk:
                tmp_keys = adata.uns[pk].keys()
                for i in tmp_keys:
                    #if type(adata.uns[pk][i]) == type(dict()): continue
                    nn_key_list.append(adata.uns[pk][i])
                break # only return connectivity and distance keys for one network

    if len(nn_key_list) == 0:
        nn_key_list = None
    return nn_key_list

def extract_NN_connectivities(adata = None, key_added = None):
    ad_guard(adata)
    
    connectivities = None
    nn_key_list = find_NN_keys(adata=adata, key_added=key_added)

    if type(nn_key_list) is type(None):
        return connectivities
    
    for nk in nn_key_list:
        if "connectivities" in nk:
            connectivities = adata.obsp[nk]
    
    return connectivities


def extract_NN_distances(adata = None, key_added = None):
    ad_guard(adata)
    
    distances = None
    nn_key_list = find_NN_keys(adata=adata, key_added=key_added)

    if type(nn_key_list) is type(None):
        return distances
    
    for nk in nn_key_list:
        if "distances" in nk:
            distances = adata.obsp[nk]
    
    return distances

def extract_NN_info(adata = None, key_added = None):
    ad_guard(adata)
    nn_keys = find_NN_keys(adata, key_added=key_added)
    nn_info = None
    for nk in nn_keys:
        if type(nk) is dict:
            nn_info = pd.Series(nk)
    return nn_info
