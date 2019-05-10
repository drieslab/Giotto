#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:09:27 2019

@author: rubendries
"""

import igraph as ig
import leidenalg as la
import pandas as pd
import networkx as nx

def python_leiden(df, partition_type, initial_membership=None, weights=None, n_iterations=2, seed=None, resolution_parameter = 1):
    
    # create networkx object
    Gx = nx.from_pandas_edgelist(df = df, source = 'from', target =  'to', edge_attr = 'weight')  
    
    # get weight attribute
    myweights = nx.get_edge_attributes(Gx, 'weight')
    
    # convert to igraph
    G = ig.Graph.TupleList(Gx.edges(), directed=False)
    G.es['weight'] = list(myweights.values())


    if partition_type == 'RBConfigurationVertexPartition':
        partition = la.find_partition(G, partition_type=la.RBConfigurationVertexPartition, initial_membership=initial_membership, weights=weights, n_iterations=n_iterations, seed=seed, resolution_parameter = resolution_parameter)
    elif partition_type == 'ModularityVertexPartition':
        partition = la.find_partition(G, partition_type=la.ModularityVertexPartition, initial_membership=initial_membership, weights=weights, n_iterations=n_iterations, seed=seed)
    else:
        print('no other configurations have been tested')

    # create dataframe with results
    
    vname = partition.graph.vs['name']
    membership = partition.membership
    membership_plus1 = [x+1 for x in membership]
    datadict = {'V' : vname, 'mem' : membership_plus1}
    leiden_dfr = pd.DataFrame(datadict)
    leiden_dfr = leiden_dfr.set_index('V')
    
    leiden_dfr = pd.DataFrame(datadict)

    return(leiden_dfr)
    
    
