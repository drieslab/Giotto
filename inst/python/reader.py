#!/usr/bin/python
from smfishHmrf.HMRFInstance import HMRFInstance
from smfishHmrf.DatasetMatrix import DatasetMatrix, DatasetMatrixSingleField, DatasetMatrixMultiField
from smfishHmrf.spatial import rank_transform_matrix, calc_silhouette_per_gene
import sys
import os
import math
import numpy as np
import scipy
import scipy.stats
from scipy.stats import zscore
from scipy.spatial.distance import euclidean, squareform, pdist
import smfishHmrf.reader as reader
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from smfishHmrf.bias_correction import calc_bias_moving, do_pca, plot_pca
from scipy.cluster.vq import kmeans2
import argparse

def read_centroid(n):
	f = open(n)
	num_cell = 0
	for l in f:
		l = l.rstrip("\n")
		num_cell+=1
	f.close()
	Xcen = np.empty((num_cell, 2), dtype="float32")
	field = np.empty((num_cell), dtype="int32")
	f = open(n)
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		x1, x2 = float(ll[0]), float(ll[1])
		t_id = int(ll[-1].split("_")[1]) - 1
		Xcen[t_id, :] = [x1, x2]
		field[t_id] = 100
	f.close()
	return Xcen, field

def read_graph(n):
	f = open(n)
	edges = set([])
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		e1, e2 = ll
		e1_id = int(e1.split("_")[1])-1
		e2_id = int(e2.split("_")[1])-1
		edges.add(tuple(sorted([e1_id, e2_id])))
	f.close()
	return edges

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="HMRF.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-l", "--location", dest="location", type=str, required=True)
	parser.add_argument("-g", "--genes", dest="genes", type=str, required=True)
	parser.add_argument("-n", "--network", dest="network", type=str, required=True)
	parser.add_argument("-e", "--expression", dest="expression", type=str, required=True)
	parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True)
	parser.add_argument("-a", "--name", dest="name", type=str, required=True)
	
	parser.add_argument("-k", "--k", dest="k", type=int, required=True)
	parser.add_argument("-b", "--betas", help="three numbers: start_beta, beta_increment, num_beta (e.g. 0 2.0 50)", nargs=3, dest="betas", type=float, required=True)
	parser.add_argument("-t", "--tolerance", dest="tolerance", type=float, help="tolerance value", default=1e-10)
	parser.add_argument("-z", "--zscore", type=str, dest="zscore", choices=["rowcol", "colrow", "none"], default="none", help="zscore the matrix after subsetting to spatial genes. Rowcol: row(gene) first, column(cell) next.")
	parser.add_argument("-i", "--numinit", type=int, dest="num_init", default=100, help="number of initializations")

	args = parser.parse_args()

	#print args
	#sys.exit(0)	
	
	Xcen, field = read_centroid(args.location)
	mat = pd.read_table(args.expression, sep=" ", header=0, index_col=0)
	
	#print mat.index	
	genes = []
	for g in range(mat.index.shape[0]):
		genes.append(str(mat.index[g]))
	#print genes
	genes_good = reader.read_genes(args.genes)
	expr = np.copy(mat.values)

	new_dset = DatasetMatrixSingleField(expr, genes, None, Xcen)	
	edges = read_graph(args.network)
	points = set([])
	adjacent = {}
	for e1,e2 in edges:
		points.add(e1)
		points.add(e2)
	ncell = expr.shape[1]
	ngene = expr.shape[0]
	#print ncell, ngene
	
	dist = pdist(Xcen, metric="euclidean")
	dist = squareform(dist)
	for i in range(ncell):
		if i in points: continue
		dist_i = sorted([(dist[i,j],j) for j in range(ncell) if i!=j])
		edges.add(tuple(sorted([i, dist_i[0][1]])))
	for e1,e2 in edges:
		adjacent.setdefault(e1, set([]))
		adjacent.setdefault(e2, set([]))
		adjacent[e1].add(e2)
		adjacent[e2].add(e1)
	new_dset.edges = edges
	new_dset.adjacent = adjacent
	new_dset.calc_independent_region()

	t_dset = new_dset.subset_genes(genes_good)

	if args.zscore=="colrow":
		t_dset.expr = zscore(t_dset.expr, axis=0) #per column (cell)
		t_dset.expr = zscore(t_dset.expr, axis=1) #per row (gene)
	elif args.zscore=="rowcol":
		t_dset.expr = zscore(t_dset.expr, axis=1) #per row (gene)
		t_dset.expr = zscore(t_dset.expr, axis=0) #per col (cell)
		
	outdir = args.outdir
	st_beta, incr_beta, num_beta = args.betas
	st_beta = float(st_beta)
	incr_beta = float(incr_beta)
	num_beta = int(num_beta)
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	this_hmrf = HMRFInstance(args.name, outdir, t_dset, args.k, st_beta, incr_beta, num_beta, tolerance=args.tolerance)
	this_hmrf.init(nstart=args.num_init)
	this_hmrf.run()
