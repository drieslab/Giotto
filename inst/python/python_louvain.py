import community
import networkx as nx
import pandas as pd

def python_louvain(df, resolution, randomize=None, random_state=None):
    G = nx.from_pandas_edgelist(df = df, source = 'from', target =  'to', edge_attr = 'weight')
    partition = community.best_partition(graph=G, resolution=resolution, weight='weight',
    randomize=randomize, random_state= random_state)
    louvain_dfr = pd.DataFrame.from_dict(data=partition, orient='index')
    return(louvain_dfr)
    


