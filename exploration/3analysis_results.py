#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import glob
import seaborn as sns
import matplotlib.pyplot as plt

import os
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)


def ranked(data1, data2) :
    dist = 0
    dis1 = list(data1['node'])
    dis2 = list(data2['node'])
    dico_rk1 = dict()
    dico_dis1 = dict()
    dico_rk2 = dict()
    dico_dis2 = dict()
    for k in range(len(data1)) :
        dico_rk1[data1['node'][k]] = k
        dico_dis1[data1['node'][k]] = data1['score'][k]
        dico_rk2[data2['node'][k]] = k
        dico_dis2[data2['node'][k]] = data2['score'][k]
    for k in range(len(dis1)) : 
        dist1 = abs(dico_rk1[dis1[k]] - dico_rk2[dis1[k]])
        dist2 = abs(dico_rk1[dis2[k]] - dico_rk2[dis2[k]])
        dist += np.sqrt(dist1**2 + dist2**2)*(1/(((dico_rk1[dis1[k]]+dico_rk1[dis2[k]]+1)/2)**2))
    return dist
    


name = 'results/param_'
num = 107
multi_1 = list()
multi_2 = list()
multi_3 = list()
for k in range(0, num) :
    path = name + str(k) + '/multiplex_protein.tsv'
    temp = pd.read_csv(path, sep = '\t')
    temp = temp[['node', 'score']]
    if (np.isnan(temp['score'][0]) == False) and (temp['score'][0] != 0) :
        multi_1.append(temp)
    path = name + str(k) + '/multiplex_disease.tsv'
    temp = pd.read_csv(path, sep = '\t')
    temp = temp[['node', 'score']]
    if (np.isnan(temp['score'][0]) == False) and (temp['score'][0] != 0) :
        multi_2.append(temp)
    path = name + str(k) + '/multiplex_drug.tsv'
    temp = pd.read_csv(path, sep = '\t')
    temp = temp[['node', 'score']]
    if (np.isnan(temp['score'][0]) == False) and (temp['score'][0] != 0) :
        multi_3.append(temp)
    

##############################################################################
## heatmap :

heat = np.zeros((len(multi_1), len(multi_1)), dtype = object)
for i in range(len(multi_1)) :
    for j in range(len(multi_1)) :
        heat[i,j] = ranked(multi_1[i], multi_1[j])
        print(i,j)
np.save('heat', heat)

heat2 = np.zeros((len(multi_2), len(multi_2)), dtype = object)
for i in range(len(multi_2)) :
    for j in range(len(multi_2)) :
        heat2[i,j] = ranked(multi_2[i], multi_2[j])
        print(i,j)
np.save('heat2', heat2)

heat3 = np.zeros((len(multi_3), len(multi_3)), dtype = object)
for i in range(len(multi_3)) :
    for j in range(len(multi_3)) :
        heat3[i,j] = ranked(multi_3[i], multi_3[j])
        print(i,j)
np.save('heat3', heat3)


heat = np.load('heat.npy', allow_pickle=True)
heat_df = pd.DataFrame(heat, dtype=float)
sns.heatmap(heat_df, cmap="Blues")
plt.savefig('heatmap.png', dpi=1000)
plt.close()

heat2 = np.load('heat2.npy', allow_pickle=True)
heat2_df = pd.DataFrame(heat2, dtype=float)
sns.heatmap(heat2_df, cmap="Blues")
plt.savefig('heatmap2.png', dpi=1000)
plt.close()

heat3 = np.load('heat3.npy', allow_pickle=True)
heat3_df = pd.DataFrame(heat3, dtype=float)
sns.heatmap(heat3_df, cmap="Blues")
plt.savefig('heatmap3.png', dpi=1000)
plt.close()
