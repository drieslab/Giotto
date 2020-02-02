#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:13:31 2019

@author: rubendries
"""

import scipy
import scipy.stats
import sys
import re
import os
import numpy as np
import math
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
#import smfishHmrf.reader as reader
import pandas as pd

def get_distance_per_FD_2(mr_dissimilarity_FD, num_cell, clust, outcome=[1,2]):
	c1 = np.where(clust==1)[0]
	c2 = np.where(clust==2)[0]
	within_dist = mr_dissimilarity_FD[np.ix_(c1, c1)]
	across_dist = mr_dissimilarity_FD[np.ix_(c1, c2)]
	mm_vec = (np.sum(within_dist, axis=1) - within_dist.diagonal()) / float(within_dist.shape[0] - 1)
	mn_vec = np.mean(across_dist, axis=1)
	sil_vec = (mn_vec - mm_vec)/np.max(np.concatenate(([mn_vec], [mm_vec])), axis=0)
	avg_clust1_sil = np.mean(sil_vec)
	return avg_clust1_sil

def rank_transform_matrix(mat, rbp_p = 0.99, reverse=True):
	dim1 = mat.shape[0]
	dim2 = mat.shape[1]
	rank_forward = np.empty([dim1, dim2])
	for c1 in range(dim1):
		rd = scipy.stats.rankdata(mat[c1,:])
		if reverse==True:
			rd = dim2 - rd + 1
		rank_forward[c1, :] = rd
	rank_backward = np.empty([dim1, dim2])
	for c1 in range(dim2):
		rd = scipy.stats.rankdata(mat[:,c1])
		if reverse==True:
			rd = dim1 - rd + 1
		rank_backward[:, c1] = rd
	mutual_rank_rbp = np.empty([dim1, dim2])
	mutual_rank = np.empty([dim1, dim2])
	for c1 in range(dim1):
		for c2 in range(dim2):
			ma = math.sqrt(rank_forward[c1, c2] * rank_backward[c1, c2])
			mutual_rank_rbp[c1, c2] = (1-rbp_p) * math.pow(rbp_p, ma - 1) 
			mutual_rank[c1, c2] = ma
	dissimilarity = np.empty([dim1, dim2])
	for c1 in range(dim1):
		for c2 in range(dim2):
			dissimilarity[c1, c2] = 1 - mutual_rank_rbp[c1, c2] / (1-rbp_p)
	return dissimilarity

def calc_silhouette_per_gene(genes=None, expr=None, dissim=None, examine_top=0.1):
	if genes is None or expr is None or dissim is None:
		sys.stderr.write("Need genes, expr, dissim\n")
		return ;
	sys.stdout.write("Started" + "\n")
	sil = []
	ncell = expr.shape[1]
	ex = int((1.0-examine_top)*100.0)
	for ig,g in enumerate(genes):
		sys.stdout.write("%s %d / %d\n" % (g, ig, len(genes)))
		cutoff = np.percentile(expr[ig,:], ex)
		clust = np.zeros((ncell), dtype="int32")
		clust[np.where(expr[ig,:]>=cutoff)[0]] = 1
		clust[np.where(expr[ig,:]<cutoff)[0]] = 2
		#avg_sil_rank, all_silhouette = get_distance_per_FD(dissim, ncell, clust)
		#avg_sil_rank, all_silhouette = get_distance_per_FD_2(dissim, ncell, clust, outcome=[1,2])
		avg_clust1_sil = get_distance_per_FD_2(dissim, ncell, clust, outcome=[1,2])
		#sil.append((g, avg_sil_rank, np.mean(all_silhouette[np.where(clust==1)[0]])))
		#sil.append((g, avg_sil_rank, avg_clust1_sil))
		sil.append((g, -1, avg_clust1_sil))
	res = []
	for ig,g in enumerate(genes):
		this_avg = sil[ig][1]
		this_sil = sil[ig][2]
		res.append((g, this_sil))
	#res.sort(lambda x,y:cmp(x[1], y[1]), reverse=True)
	res.sort(key=itemgetter(1), reverse=True)
	return res


def python_spatial_genes(spatial_locations, expression_matrix,
                         metric = "euclidean",
                         rbp_p = 0.95, examine_top = 0.3):
    
    Xcen =  spatial_locations
    mat = expression_matrix
    genes = []
    
    for g in range(mat.index.shape[0]):
        genes.append(str(mat.index[g]))
    expr = np.copy(mat.values)
    
    ncell = Xcen.shape[0] 
    sys.stdout.write("Calculate all pairwise Euclidean distance between cells using their physical coordinates\n")
    euc = squareform(pdist(Xcen, metric=metric))
    sys.stdout.write("Rank transform euclidean distance, and then apply exponential transform\n")
    dissim = rank_transform_matrix(euc, reverse=False, rbp_p=rbp_p)
    sys.stdout.write("Compute silhouette metric per gene\n")
    res = calc_silhouette_per_gene(genes=genes, expr=expr, dissim=dissim, examine_top=0.3)
    
    return res
    
