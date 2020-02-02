#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:13:31 2019

@author: Qian Zhu
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
	
	print("Start ranking forward...")
	for c1 in range(dim1):
		rd = scipy.stats.rankdata(mat[c1,:])
		if reverse==True:
			rd = dim2 - rd + 1
		rank_forward[c1, :] = rd
		if c1%1000==0:
			print("Done %d" % c1)
	
	#rank_forward = scipy.stats.mstats.rankdata(mat, axis=0)
	print("Finished ranking forward...")
	'''
	if reverse==True:
		print("Start adjusting forward...")
		rank_forward = np.add(np.subtract(dim2, rank_forward), 1)
		print("Finished adjusting forward...")
	'''
	
	rank_backward = np.empty([dim1, dim2])
	
	print("Start ranking backward...")
	for c1 in range(dim2):
		rd = scipy.stats.rankdata(mat[:,c1])
		if reverse==True:
			rd = dim1 - rd + 1
		rank_backward[:, c1] = rd
		if c1%1000==0:
			print("Done %d" % c1)
	
	#rank_backward = scipy.stats.mstats.rankdata(mat, axis=1)
	print("Finished ranking backward...")
	'''
	if reverse==True:
		print("Start adjusting backward...")
		rank_backward = np.add(np.subtract(dim1, rank_backward), 1)
		print("Finished adjusting backward...")
	'''
	mutual_rank_rbp = np.empty([dim1, dim2])
	mutual_rank = np.empty([dim1, dim2])

	
	print("Calculate mutual rank...")
	ma = np.sqrt(np.multiply(rank_forward, rank_backward))
	print("Calculate exponential transform...")
	mutual_rank_rbp = np.multiply(1-rbp_p, np.power(rbp_p, np.subtract(ma, 1)))
	print("Finished exponential transform...")
	mutual_rank = ma
	
	
	dissimilarity = np.empty([dim1, dim2])
	'''
	print("Calculate dissimilarity...")
	for c1 in range(dim1):
		for c2 in range(dim2):
			ma = math.sqrt(rank_forward[c1, c2] * rank_backward[c1, c2])
			mutual_rank_rbp[c1, c2] = (1-rbp_p) * math.pow(rbp_p, ma - 1) 
			mutual_rank[c1, c2] = ma
			dissimilarity[c1, c2] = 1 - mutual_rank_rbp[c1, c2] / (1-rbp_p)
		if c1%1000==0:
			print("Done %d" % c1)
	'''
	'''
	for c1 in range(dim1):
		for c2 in range(dim2):
			dissimilarity[c1, c2] = 1 - mutual_rank_rbp[c1, c2] / (1-rbp_p)
	'''
	print("Calculate dissimilarity...")
	dissimilarity = np.subtract(1, np.divide(mutual_rank_rbp, 1-rbp_p))
	print("Finished dissimilarity...")
	return dissimilarity

def calc_silhouette_per_gene(genes=None, expr=None, dissim=None, examine_top=0.1, seed=-1):
	if genes is None or expr is None or dissim is None:
		sys.stderr.write("Need genes, expr, dissim\n")
		return ;
	if seed!=-1 and seed>=0:
		np.random.seed(seed)
	sys.stdout.write("Started 2 " + "\n")
	sil = []
	ncell = expr.shape[1]
	ex = int((1.0-examine_top)*100.0)
	for ig,g in enumerate(genes):
		cutoff = np.percentile(expr[ig,:], ex)
		clust = np.zeros((ncell), dtype="int32")
		gt_eq = np.where(expr[ig,:]>=cutoff)[0]
		lt = np.where(expr[ig,:]<cutoff)[0]
		if gt_eq.shape[0]>int(ncell*examine_top):
			num_filter = gt_eq.shape[0] - int(ncell*examine_top)
			ss = np.random.choice(gt_eq, size=num_filter, replace=False)
			clust[gt_eq] = 1
			clust[lt] = 2
			clust[ss] = 2
		elif gt_eq.shape[0]<int(ncell*examine_top):
			num_filter = int(ncell*examine_top) - gt_eq.shape[0]
			ss = np.random.choice(lt, size=num_filter, replace=False)
			clust[gt_eq] = 1
			clust[lt] = 2
			clust[ss] = 1
		else:
			clust[gt_eq] = 1
			clust[lt] = 2
		'''
		if cutoff==0:
			val_gt = np.where(expr[ig,:]>0)[0]
			val_iszero = np.where(expr[ig,:]==0)[0]
			num_filler = int(ncell*examine_top) - val_gt.shape[0]
			ss = np.random.choice(val_iszero, size=num_filler, replace=False)
			clust[val_iszero] = 2
			clust[val_gt] = 1
			clust[ss] = 1
		'''
		'''
		'''	
		sys.stdout.write("%s %d / %d\n" % (g, ig, len(genes)))
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
    res = calc_silhouette_per_gene(genes=genes, expr=expr, dissim=dissim, examine_top=examine_top)
    
    return res
    
