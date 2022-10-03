#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
from validation import protocol1, protocol2, loocv
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
import glob 

# define set of nodes of the interesting multiplex :
multiplex2_path = glob.glob('networks/multiplex/2/' + '*.gr')
multiplex2_nodes = set()
for k in range(len(multiplex2_path)) :
    temp = pd.read_csv(multiplex2_path[k], sep ='\t', header = None)
    temp = temp.drop(temp[temp[0] == temp[1]].index).reset_index()
    temp_net = nx.from_pandas_edgelist(temp, 0, 1)
    multiplex2_nodes = multiplex2_nodes.union(set(temp_net.nodes))

# define set of nodes of the interacted multiplex
multiplex1_path = glob.glob('networks/multiplex/1/' + '*.gr')
multiplex1_nodes = set()
for k in range(len(multiplex1_path)) :
    temp = pd.read_csv(multiplex1_path[k], sep ='\t', header = None)
    temp = temp.drop(temp[temp[0] == temp[1]].index).reset_index()
    temp_net = nx.from_pandas_edgelist(temp, 0, 1)
    multiplex1_nodes = multiplex1_nodes.union(set(temp_net.nodes))
    
# choose the bipartite :
bip = '1_2'
bipartite = pd.read_csv(bip + '.gr', sep ='\t', header = None)
bipartite = bipartite.drop_duplicates()

# remove nodes only in bipartite and note in multiplexes
if (bip == '1_2') :
    bipartite = bipartite[bipartite[0].isin(multiplex1_nodes)]
    bipartite = bipartite[bipartite[1].isin(multiplex2_nodes)]
    multiplex2_nodes = multiplex2_nodes.intersection(set(bipartite[1]))
else :  # bip == '2_1'
    bipartite = bipartite[bipartite[1].isin(multiplex1_nodes)]
    bipartite = bipartite[bipartite[0].isin(multiplex2_nodes)]
    multiplex2_nodes = multiplex2_nodes.intersection(set(bipartite[0]))

# process of evaluation :
# choose the method of evaluation; loocv, protocol1 (modified \
# link prediction), protocol2 (link prediction)
num_cpu = 22
results = loocv(path, multiplex2_nodes, bipartite, bip, num_cpu)
results.to_csv("evaluation.tsv", sep = '\t', header = None, index = False)              
        
plt.figure()
plt.plot(ECDF(list(results['rank'])).x, ECDF(list(results['rank'])).y, 'red', label = "multi-3")
plt.legend()
plt.savefig('rank.png', format = 'png', dpi = 1000)
        
