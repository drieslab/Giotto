import sys
import os
import re
import numpy as np
import subprocess
import math
import scipy
import silhouetteRank.spatial_genes as spatial_genes
from shutil import copyfile
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
from sklearn.metrics import roc_auc_score
import argparse
import silhouetteRank
import silhouetteRank.prep as prep
import silhouetteRank.evaluate_exact_one_2b as evaluate_exact_one_2b
import silhouetteRank.use_previous_cluster as use_previous_cluster
import silhouetteRank.combine as combine

def silhouette_rank(expr="expression.txt", centroid="Xcen.good", overwrite_input_bin=False,
rbp_ps=[0.95, 0.99], examine_tops=[0.005, 0.010, 0.050, 0.100, 0.300], matrix_type="dissim",
logdir="./logs", num_core=4, parallel_path="/usr/bin", output=".", query_sizes=10):
	args = argparse.Namespace(expr=expr, centroid=centroid, rbp_ps=rbp_ps, examine_tops=examine_tops,
	matrix_type=matrix_type, output=output, query_sizes=query_sizes, overwrite_input_bin=overwrite_input_bin,
	logdir=logdir, parallel_path=parallel_path, num_core=num_core)

	if not os.path.isdir(args.output):
		os.mkdir(args.output)
	if not os.path.isdir(args.logdir):
		os.mkdir(args.logdir)

	
	args1 = argparse.Namespace(expr=args.expr, centroid=args.centroid, rbp_ps=args.rbp_ps, examine_tops=args.examine_tops, matrix_type=args.matrix_type, output=args.output, query_sizes=args.query_sizes, overwrite_input_bin=args.overwrite_input_bin)
	prep.do_one(args1)

	fw = open("%s/args" % args.output, "w")
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			freq_file = "%s/result_5000_%.2f_%.3f/gene.freq.good.txt" % (args.output, rbp_p, examine_top)
			if args.matrix_type=="sim":
				freq_file = "%s/result_sim_5000_%.2f_%.3f/gene.freq.good.txt" % (args.output, rbp_p, examine_top)
			uniq_freq = 0
			f = open(freq_file)
			for l in f:
				l = l.rstrip("\n")
				uniq_freq+=1
			f.close()
			num_query_sizes = args.query_sizes
			if uniq_freq<=num_query_sizes:
				num_query_sizes = uniq_freq
			for i in range(num_query_sizes):
				fw.write("%.2f\n" % rbp_p)
				fw.write("%.3f\n" % examine_top)
				fw.write("%d\n" % i)
	fw.close()

	fw = open("%s/args.basic" % args.output, "w")
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			fw.write("%.2f\n" % rbp_p)
			fw.write("%.3f\n" % examine_top)
	fw.close()

	bin_path = os.path.dirname(silhouetteRank.__file__)
	for i in range(4):
		bin_path = os.path.dirname(bin_path)
	bin_path = os.path.join(bin_path, "bin")	

	sys.stderr.write("Start calculating silhouette rank...\n")
	sys.stderr.flush()
	cmd = "cat '%s'/args.basic | '%s'/parallel --jobs %d --max-args=2 \\''%s'\\'''/silhouette_rank_main -x \\''%s'\\''' -c \\''%s'\\''' -r {1} -e {2} -m %s -o \\''%s'\\''' \"2>\" \\''%s'\\'''/real_{1}_{2}.out" % (args.output, args.parallel_path, args.num_core, bin_path, args.expr, args.centroid, args.matrix_type, args.output, args.logdir)
	os.system(cmd)

	sys.stderr.write("Start random...\n")
	sys.stderr.flush()
	cmd="cat '%s'/args | '%s'/parallel --jobs %d --max-args=3 \\''%s'\\'''/silhouette_rank_random -r {1} -e {2} -m %s -o \\''%s'\\''' -q {3} \"2>\" \\''%s'\\'''/{1}_{2}_{3}.out" % (args.output, args.parallel_path, args.num_core, bin_path, args.matrix_type, args.output, args.logdir)
	os.system(cmd)

	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			random_dir = "%s/result_sim_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
			output_score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.output, rbp_p, examine_top)
			if args.matrix_type=="dissim":
				random_dir = "%s/result_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
				score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
				output_score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.output, rbp_p, examine_top)
			args1 = argparse.Namespace(expr=args.expr, centroid=args.centroid, examine_top=examine_top, input=score_file, input_random=random_dir, output=output_score_file, outdir=args.output, query_sizes=args.query_sizes, overwrite_input_bin=args.overwrite_input_bin)
			use_previous_cluster.do_one(args1)

	combined_file = "%s/silhouette.overall.pval.txt" % args.output
	if args.matrix_type=="sim":
		combined_file = "%s/silhouette.sim.overall.pval.txt" % args.output
	args1 = argparse.Namespace(rbp_ps=args.rbp_ps, examine_tops=args.examine_tops, matrix_type=args.matrix_type, input=args.output, output=combined_file)
	combine.do_one(args1)
	 
	res = []
	f = open(combined_file)
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		res.append((ll[0], float(ll[1]), float(ll[2])))
	f.close()
	return res
