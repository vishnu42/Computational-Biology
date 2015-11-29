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

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
from scipy.spatial.distance import *
from pyemd import emd
from os import listdir
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

def read_graphs(indir, output, k):

	index = 0
	#Reference Dataset
	fp = open(output,'w')

	dic_list = []
	ground_data = []
	test_data = []
	#print dataset

	folder_name = indir + '/' 

	filenames = [f for f in listdir(folder_name)]
	for filename_out in filenames:
		model_name_o = filename_out.split('_')[0]
		Graph_List.append(folder_name+filename_out)
		G_1 = nx.read_edgelist(folder_name+filename_out)
    
		for filename_inner in filenames:
			key_1 = filename_out + filename_inner
			key_2 = filename_inner + filename_out
			if not(key_1 in dic_list) and not(key_2 in dic_list):
				dic_list.append(key_1)
				dic_list.append(key_2)
				model_name_i = filename_inner.split('_')[0]
				G_2 = nx.read_edgelist(folder_name+filename_inner)
	   
           			if model_name_i == model_name_o:
                			Reference_matrix[filename_out+':'+filename_inner] = 1

		                        Our_matrix[filename_out+':'+filename_inner] = 1
	        			ground_data.append(1)
					if not(filename_inner==filename_out):

						fp.write(filename_out)
	        				fp.write(filename_inner)
	        				v = com_sig_dist(G_1, G_2, k)
	        				test_data.append(v)
	        				print 'match',str(v)
						fp.write('match'+str(v)+'\n')
				else:
					fp.write(filename_out)
				        fp.write(filename_inner)
                			Reference_matrix[filename_out+':'+filename_inner] = 0
	       				v = com_sig_dist(G_1, G_2, k)
					Our_matrix[filename_out+':'+filename_inner] = v
					ground_data.append(0)
					test_data.append(v)
					print 'mismatch',str(v)
					fp.write('mismatch'+str(v)+'\n')
	        

	ground_data = np.array(ground_data)
	test_data = np.array(test_data)
	
	'''
	prec,rec,_ = precision_recall_curve(ground_data,ground_score)
	plt.clf()
	plt.plot(rec,prec, label='Prec-recall curve')
	plt.savefig('c.png')
	
	for x in Our_matrix.keys():
		print x, Our_matrix[x]
	'''	
	fp.close()

    
    
#===========================================================================#
#                 Compute signature for G                                   #
#===========================================================================#

    
def covariance_sig(G, k):   

	adj_matrix = nx.adjacency_matrix(G)
	num_nodes = G.number_of_nodes();

	e1 = np.zeros((num_nodes, 1))+1

	power_iter_matrix = (G.number_of_nodes() * (adj_matrix * e1))/np.linalg.norm(adj_matrix * e1)
	for i in range(2,k+1):
		p_iter = np.linalg.matrix_power(adj_matrix, i) * e1
		p_iter = (G.number_of_nodes() * (p_iter))/np.linalg.norm(p_iter)
		power_iter_matrix = np.column_stack([power_iter_matrix, p_iter])
	
	cov_matrix = [[0 for x in range(k)] for x in range(k)]
	for i in range(num_nodes):
		C_i = np.matrix.transpose(power_iter_matrix[i,0:]-1) * (power_iter_matrix[i,0:]-1)
		cov_matrix += C_i
	
	return cov_matrix
	
#===========================================================================#
#                 Distance between G1,G2                                    #
#===========================================================================#
    
def com_sig_dist(G_1, G_2, k):
	cov_matrix_A = covariance_sig(G_1, k)
	cov_matrix_B = covariance_sig(G_2, k)
	
	det_AM_of_cov_matrices = np.linalg.det((cov_matrix_A+cov_matrix_B)/2)
	sqrt_det_covA_covB = math.sqrt(np.linalg.det(cov_matrix_A)*np.linalg.det(cov_matrix_B))

	dist = np.log(det_AM_of_cov_matrices/sqrt_det_covA_covB)/2
	sim = np.exp(-dist)

	return dist

def main():
	parser = argparse.ArgumentParser(description="Compute comsig distances between graphs")
	parser.add_argument('--input', help="input directory")
	parser.add_argument('--output', help="output file")
	parser.add_argument('--k', type=int, default=5, help="number of iterations")
	#parser.add_argument('--threads', type=int, default=10, help="number of threads to use")
	args = parser.parse_args()
	#read_graphs(args.input, args.output, args.k, args.threads)
	read_graphs(args.input, args.output, args.k)

if __name__ == "__main__":
	main()
