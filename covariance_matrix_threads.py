from __future__ import division
from optparse import OptionParser
from collections import deque
import math
from math import log
import argparse
import os
import sys,time,logging,random
import snap
import numpy as np
import matplotlib
from subprocess import call
import networkx as nx
import csv
import glob
from multiprocessing import Pool

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
from scipy.spatial.distance import *
from pyemd import emd
from os import listdir
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve


#from igraph import *
matplotlib.use('Agg')
#from graphviz import *
Graph_List = []
models = ['ER','GEO','PAM','DDM']
rho = [0.005,0.01]
nodes = [1000]
Reference_matrix = {}
Our_matrix = {}
#dataset = []
k_value = 0

def read_graphs(indir, output, threads):

	index = 0
	#Reference Dataset
	fp = open(output,'w')

	p = Pool(processes=threads)
	paths = glob.glob(indir + "/*.txt")
	#print paths
	sig_dic = {}
	
	for res in p.imap_unordered(covariance_sig, paths, chunksize=1):
        	print("processed {}".format(res[0]))
        	sig_dic[res[0]] = [res[1],res[2]]
		#print res[0],sig_dic[res[0]]

	p.close()
	p.join()

	dic_list = []
	ground_data = []
	test_data = []
	#print dataset
	flag_set = set([])
	validTuples = []
	for k1, v1 in sig_dic.iteritems():
		for k2, v2 in sig_dic.iteritems():
		#start computing if densities are equal
			#if True:#k1.Density == k2.Density:
			if k1 != k2: 
				validTuples.append((k1, k2))
	#print "Valid tuples: ", validTuples
	fp.write("Network1,\t\t\t Network2,\t\t\t dist, Sim, G-sim\n")
	p = Pool(processes=threads)
	for r in p.imap_unordered(com_sig_dist, ((t[0], t[1], sig_dic[t[0]][0], sig_dic[t[1]][0], sig_dic[t[0]][1], sig_dic[t[1]][1]) for t in validTuples), chunksize=10):
        	k1, k2 = r[0], r[1]
        	dis, sim, g_sim  = r[2], r[3], r[4]
		ground_data.append(g_sim)
		test_data.append(sim)
        	print("Done G1: {}, G2: {}".format(k1, k2))
        	fp.write("{}, {}, {:.5f}, {:.5f}, {}\n".format(k1, k2, dis, sim, g_sim))
	
	p.close()
	p.join()
	fp.close()

	ground_data = np.array(ground_data)
	test_data = np.array(test_data)

	'''for i in range(2,10):
		new_data = []
		for sim in test_data:
			if (sim > float(i/10)):
				sim = 1
				new_data.append(1)
			else:
				sim = 0
				new_data.append(0)
	'''
	prec,rec,_ = precision_recall_curve(ground_data,test_data)
	print "precision =", prec
	print "recall =", rec
	plt.clf()
	plt.plot(rec,prec, label='Prec-recall curve')
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.savefig('covariance_prec-recall.png')
	AUPR = average_precision_score(ground_data, test_data)
	print "Area under prec-recall curve: ", AUPR
	#plt.plot(prec,rec, label='recall-prec curve')
	#plt.savefig('cov_'+str(i)+'.png')
	#plt.savefig('covariance_recall_prec.png')

def mygraph():
	
	ground_data = [1,1,0,1,0,0,1,0]
	test_data = [0,1,0,0,0,0,0,1]

	print "++++ ground_data: ", ground_data
	print "++++ test_data: ", test_data	
	prec,rec,_ = precision_recall_curve(ground_data,test_data)
	plt.clf()
	plt.plot(rec,prec, label='Prec-recall curve')
	plt.savefig('my.png')
 
#===========================================================================#
#                 Compute signature for G                                   #
#===========================================================================#

def covariance_sig(fn):   
	#print "k value", k_value
	#print "In cov sig", fn
	filename = fn.split('/')[-1]
	#print 'filename =', filename
	model_name = filename.split('_')[0]
	#print 'modelname =', model_name
	G = nx.read_edgelist(fn)
	adj_matrix = nx.adjacency_matrix(G)
	num_nodes = G.number_of_nodes();
	#print "In cov sig", num_nodes

	e1 = np.zeros((num_nodes, 1))+1

	p_iter = adj_matrix * e1
	power_iter_matrix = (G.number_of_nodes() * (p_iter))/np.linalg.norm(p_iter)
	for i in range(2,k_value+1):
		p_iter = adj_matrix * p_iter
		p_iter = (G.number_of_nodes() * (p_iter))/np.linalg.norm(p_iter)
		power_iter_matrix = np.column_stack([power_iter_matrix, p_iter])
	#print "power matrix", power_iter_matrix
	
	cov_matrix = [[0 for x in range(k_value)] for x in range(k_value)]
	for i in range(num_nodes):
		C_i = np.matrix.transpose(power_iter_matrix[i,0:]-1) * (power_iter_matrix[i,0:]-1)
		cov_matrix += C_i
	
	return (fn, cov_matrix/num_nodes, model_name)
	
#===========================================================================#
#                 Distance between G1,G2                                    #
#===========================================================================#
    
def com_sig_dist(tupArg):
	cov_matrix_A = tupArg[2] 
	#print "++ Covariance Matrix_A: ", cov_matrix_A
	cov_matrix_B = tupArg[3] 
	#print "++ Covariance Matrix_B: ", cov_matrix_B
	model_name_A = tupArg[4]
	model_name_B = tupArg[5]
	ground_sim = 0
	if (model_name_A == model_name_B):
		ground_sim = 1
	
	det_AM_of_cov_matrices = np.linalg.det((cov_matrix_A+cov_matrix_B)/2)
	#print '++ Det: ', np.linalg.det((cov_matrix_A+cov_matrix_B))

	sqrt_det_covA_covB = math.sqrt(np.linalg.det(cov_matrix_A)*np.linalg.det(cov_matrix_B))
	#print '++ Sqrt: ', sqrt_det_covA_covB 

	dist = np.log(det_AM_of_cov_matrices/sqrt_det_covA_covB)/2
	sim = np.exp(-dist)
	#print "sim = ", sim

	#print '++++++++++++++++ Dist: ', dist
	return (tupArg[0], tupArg[1], dist, sim, ground_sim)

def main():
	parser = argparse.ArgumentParser(description="Compute comsig distances between graphs")
	parser.add_argument('--input', help="input directory")
	parser.add_argument('--output', help="output file")
	parser.add_argument('--k', type=int, default=5, help="number of iterations")
	parser.add_argument('--threads', type=int, default=10, help="number of threads to use")
	args = parser.parse_args()
	global k_value 
	k_value = args.k
	#print "kvalue_main", k_value
	#read_graphs(args.input, args.output, args.k, args.threads)
	#mygraph()
	#return
	read_graphs(args.input, args.output, args.threads)

if __name__ == "__main__":
	main()
