from __future__ import division
from optparse import OptionParser
from collections import deque
from math import log
import argparse
import os
import sys,time,logging,random
import snap
#from igraph import *
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




def duplication_divergence_model(Nodes,Edges,sigma):
	edge_list = []
	node_set = set(range(0,Nodes))
	G = nx.Graph()
	G.add_node(0)
	node_set.remove(0)

	while len(G.nodes()) < Nodes:
		u = random.randint(0,len(G.nodes()))
		v = list(node_set)[0]
		node_set.remove(v)
		G.add_node(v)
		if (len(G.nodes()) == 2):
			G.add_edge(u,v)
		else:
			for w in G.neighbors(u):
				if random.random() < sigma:
					G.add_edge(v,w)
			if (len(G.neighbors(v)) == 0):
				G.add_edge(u,v)
		

	return G


def main(args):
	#model_name = args.model
	Nodes = args.N
	rho = args.rho
	Edges = int(rho * Nodes * (Nodes -1) / 2)
	seed = args.seed
	seq = args.num
	max_rewire = args.rewire
	folder = './'+'dataset_'+str(seq)+'/'

	folder_rewire = './'+'dataset_rewire_'+str(seq)
	os.system("mkdir "+folder)
	density = rho
	os.system("mkdir "+folder_rewire)
	avg_deg = int(2*Edges / Nodes)
	print avg_deg

	'''


	if model_name=='PAM':
		for i in range(1,seed+1):
			file_name = folder+'PAM_'+str(Nodes)+'_'+str(Edges)+'_'+str(i)+'.txt'
			G = nx.barabasi_albert_graph(Nodes,avg_deg,seed=i)
			nx.write_edgelist(G,file_name)



	if model_name=='ER':
		for i in range(1,seed+1):
			file_name = folder+'ER_'+str(Nodes)+'_'+str(Edges)+'_'+str(i)+'.txt'
			G = nx.erdos_renyi_graph(Nodes, density,seed=i, directed=False )
			nx.write_edgelist(G,file_name)

	if model_name=='GEO':

		for i in range(1,seed+1):
			file_name = folder+'GEO_'+str(Nodes)+'_'+str(Edges)+'_'+str(i)+'.txt'
			G = nx.random_geometric_graph(Nodes, seed=i)
			nx.write_edgelist(G,file_name)
			
	if model_name == 'DDM':
		for i in range(1,seed+1):
			file_name = folder+'DDM_'+str(Nodes)+'_'+str(Edges)+'_'+str(i)+'.txt'
			sigma = 0.5			
			G = duplication_divergence_model(Nodes,Edges,sigma)
			nx.write_edgelist(G,file_name)
	
	'''

	for i in xrange(10,max_rewire,10):
		folder_tmp = folder_rewire + '/'+'dataset_rewirerate_'+str(i)
		os.system('mkdir '+folder_tmp)


	for i in range(1,seed+1):
		file_name = folder+'PAM_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
		G = nx.barabasi_albert_graph(Nodes,avg_deg,seed=i)
		for j in xrange(10,max_rewire,10):
			num_rewire = int(j*Edges/100)
			print num_rewire
			G_tmp = nx.double_edge_swap(G,num_rewire,Edges)	
			folder_tmp = folder_rewire + '/'+'dataset_rewirerate_'+str(j)
			file_name_rewire = folder_tmp +'/'+'R_PAM_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
			nx.write_edgelist(G_tmp,file_name_rewire)
		nx.write_edgelist(G,file_name)

	
	
	for i in range(1,seed+1):
		file_name = folder+'ER_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
		G = nx.erdos_renyi_graph(Nodes, density,seed=i, directed=False )
		
		for j in xrange(10,max_rewire,10):
			num_rewire = int(j*Edges/100)
			print num_rewire
			G_tmp = nx.double_edge_swap(G,num_rewire,Edges)	
			folder_tmp = folder_rewire + '/'+'dataset_rewirerate_'+str(j)
			file_name_rewire = folder_tmp  +'/'+'R_ER_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
			nx.write_edgelist(G_tmp,file_name_rewire)
		nx.write_edgelist(G,file_name)



	for i in range(1,seed+1):
		file_name = folder+'GEO_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
		G = nx.random_geometric_graph(Nodes, 0.042)
		print nx.density(G)
		nx.write_edgelist(G,file_name)
		for j in xrange(10,max_rewire,10):
			num_rewire = int(j*Edges/100)
			print num_rewire
			G_tmp = nx.double_edge_swap(G,num_rewire,Edges)	
			folder_tmp = folder_rewire + '/'+'dataset_rewirerate_'+str(j)
			file_name_rewire = folder_tmp +'/'+'R_GEO_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
			nx.write_edgelist(G_tmp,file_name_rewire)


	for i in range(1,seed+1):
		file_name = folder+'DDM_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
		sigma = 0.5			
		G = duplication_divergence_model(Nodes,Edges,sigma)
		for j in xrange(10,max_rewire,10):
			num_rewire = int(j*Edges/100)
			print num_rewire
			G_tmp = nx.double_edge_swap(G,num_rewire,Edges)	
			folder_tmp = folder_rewire + '/'+'dataset_rewirerate_'+str(j)
			file_name_rewire = folder_tmp +'/'+'R_DDM_'+str(Nodes)+'_'+str(rho)+'_'+str(i)+'.txt'
			nx.write_edgelist(G_tmp,file_name_rewire)
		nx.write_edgelist(G,file_name)








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Snap starting')
    parser.add_argument('--N', type=int, default=100, help="Number of nodes")
    parser.add_argument('--rho', type=float, default=0.01, help="Edge density")
    parser.add_argument('--num', type=int, default=100, help="dataset number")	
    parser.add_argument('--seed', type=int, default=10, help="Number of seeds")
    parser.add_argument('--rewire', type=int, default=90, help="rewire")
    args = parser.parse_args()
    main(args)
 
