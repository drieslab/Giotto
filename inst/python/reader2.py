#!/usr/bin/python
from smfishHmrf.HMRFInstance import HMRFInstance
from smfishHmrf.DatasetMatrix import DatasetMatrix, DatasetMatrixSingleField, DatasetMatrixMultiField
from smfishHmrf.spatial import rank_transform_matrix, calc_silhouette_per_gene
import sys
import os
import math
import subprocess
import numpy as np
import scipy
import scipy.stats
from scipy.stats import zscore
from scipy.spatial.distance import euclidean, squareform, pdist
import smfishHmrf.reader as reader
import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns
from smfishHmrf.bias_correction import calc_bias_moving, do_pca, plot_pca
from scipy.cluster.vq import kmeans2
import argparse

def read_centroid(n, cells):
	map_cell = {}
	for ind,val in enumerate(cells):
		map_cell[val] = ind

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
		t_id = map_cell[ll[-1]]
		#t_id = int(ll[-1].split("_")[1]) - 1
		Xcen[t_id, :] = [x1, x2]
		field[t_id] = 100
	f.close()
	return Xcen, field

def read_graph(n, cells):
	map_cell = {}
	for ind,val in enumerate(cells):
		map_cell[val] = ind

	f = open(n)
	edges = set([])
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		e1, e2 = ll
		e1_id = map_cell[e1]
		e2_id = map_cell[e2]
		#e1_id = int(e1.split("_")[1])-1
		#e2_id = int(e2.split("_")[1])-1
		edges.add(tuple(sorted([e1_id, e2_id])))
	f.close()
	return edges

def read_expression_classic(n):
	f = open(n)
	#h = f.readline().rstrip("\n").split()
	h = f.readline().rstrip("\n")
	#header begins with space
	if h.startswith(" "):
		h = h.split()
	else: #header startswith a gene name
		h = h.split()[1:]
	num_cell = len(h)
	num_gene = 0
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		#gene = ll[0]
		num_gene+=1
	f.close()
	mat = np.empty((num_gene, num_cell), dtype="float32")
	genes = []
	cells = h
	f = open(n)
	f.readline()
	gid = 0
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		genes.append(ll[0])
		mat[gid, :] = [float(v) for v in ll[1:]]
		gid+=1
	f.close()
	return mat, genes, cells

def connected_components(edges, adjacent, points):
	visited = {}
	chains = []
	for p in sorted(list(points)):
		visited[p] = False
	for p in sorted(list(points)):
		if visited[p]==False:
			new_chain = []
			visited, new_chain = DFS(p, adjacent, visited, new_chain)
			chains.append(new_chain)
	return chains

def DFS(p, adjacent, visited, new_chain):
	visited[p] = True
	new_chain.append(p)
	for nei in sorted(list(adjacent[p])):
		if visited[nei]==False:
			visited, new_chain = DFS(nei, adjacent, visited, new_chain)
	return visited, new_chain

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
	parser.add_argument("-s", "--seed", dest="seed", type=float, help="seed for random initialization of HMRF. -1 will not fix it.", default=-1)
	parser.add_argument("-z", "--zscore", type=str, dest="zscore", choices=["rowcol", "colrow", "none"], default="none", help="zscore the matrix after subsetting to spatial genes. Rowcol: row(gene) first, column(cell) next.")
	parser.add_argument("-i", "--numinit", type=int, dest="num_init", default=100, help="number of initializations")

	args = parser.parse_args()

	sys.setrecursionlimit(50000)
	#print args
	#sys.exit(0)

	mat, genes, cells = read_expression_classic(args.expression)
	print("Done reading expression")
	Xcen, field = read_centroid(args.location, cells)
	print("Done reading location")

	#mat = pd.read_table(args.expression, sep=" ", header=0, index_col=0)

	#print mat.index
	'''
	genes = []
	for g in range(mat.index.shape[0]):
		genes.append(str(mat.index[g]))
	#print genes
	expr = np.copy(mat.values)
	'''
	genes_good = reader.read_genes(args.genes)
	expr = mat


	new_dset = DatasetMatrixSingleField(expr, genes, None, Xcen)
	edges = read_graph(args.network, cells)
	print("Done reading graph")
	points = set([])
	adjacent = {}
	for e1,e2 in edges:
		points.add(e1)
		points.add(e2)
	ncell = expr.shape[1]
	ngene = expr.shape[0]
	#print ncell, ngene

	'''
	dist = pdist(Xcen, metric="euclidean")
	dist = squareform(dist)
	for i in range(ncell):
		if i in points: continue
		dist_i = sorted([(dist[i,j],j) for j in range(ncell) if i!=j])
		edges.add(tuple(sorted([i, dist_i[0][1]])))
	'''
	for e1,e2 in edges:
		adjacent.setdefault(e1, set([]))
		adjacent.setdefault(e2, set([]))
		adjacent[e1].add(e2)
		adjacent[e2].add(e1)
	new_dset.edges = edges
	new_dset.adjacent = adjacent

	print("Start calculating independent regions")
	conn = connected_components(edges, adjacent, points)

	blocks = {}
	for ind_con,con in enumerate(conn):
		all_vert = con
		set_all_vert = set(all_vert)
		map_vert = {}
		for ind,val in enumerate(all_vert):
			map_vert[val] = ind
		print("Edges for component", ind_con)
		outdir = args.outdir
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		edge_file = os.path.join(outdir, "edges.txt")
		block_file = os.path.join(outdir, "blocks.txt")
		#edge_file = "/tmp/edges.txt"
		#block_file = "/tmp/blocks.txt"
		fw = open(edge_file, "w")
		for e1, e2 in edges:
			if e1 in set_all_vert and e2 in set_all_vert:
				fw.write("%d %d\n" % (map_vert[e1]+1, map_vert[e2]+1))
		fw.close()
		import smfishHmrf
		this_path = os.path.dirname(smfishHmrf.__file__) + "/graphColoring"
		subprocess.call("java -cp '%s' -Xmx32g -Xms32g GraphColoring '%s' '%s' '%d'" % (this_path, edge_file, block_file, args.seed), shell=True)

		f = open(block_file)
		b_ind = 0
		for l in f:
			l = l.rstrip("\n")
			ll = l.split()
			#self.blocks.append(int(ll[1]))
			blocks[all_vert[b_ind]] = int(ll[1])
			b_ind+=1
		f.close()
		#self.blocks = np.array(self.blocks)

	new_blocks = []
	for b in range(0, len(blocks.keys())):
		new_blocks.append(blocks[b])
	new_dset.blocks = np.array(new_blocks)
	print("Finished calculating independent regions")


	'''
	print("Start calculating independent region")
	new_dset.calc_independent_region()
	print("Finished calculating independent region")
	'''

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
	this_hmrf.init(nstart=args.num_init, seed=args.seed)
	this_hmrf.run()
