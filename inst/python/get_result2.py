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
	parser.add_argument("-a", "--name", dest="name", type=str, required=True)
	parser.add_argument("-r", "--result", dest="result_dir", type=str, required=True)
	parser.add_argument("-k", "--k", dest="k", help="which K's result to output", type=int, required=True)
	parser.add_argument("-b", "--beta", help="which beta's result to output", dest="beta", type=float, required=True)

	args = parser.parse_args()

	cl = reader.read_classes("%s/k_%d/f%s.beta.%.1f.prob.txt" % (args.result_dir, args.k, args.name, args.beta))
	cl = cl[1:]
	for i in range(cl.shape[0]):
		sys.stdout.write("%d\n" % cl[i])
