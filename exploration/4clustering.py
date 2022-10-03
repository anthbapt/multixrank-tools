#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from collections import Counter
import itertools
import pandas as pd
import numpy as np
import math
import glob

import os
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)

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
        
    
len1 = len(multi_1[0])
len2 = len(multi_2[0])
len3 = len(multi_3[0])
heat = np.load('heat.npy', allow_pickle=True)
heat2 = np.load('heat2.npy', allow_pickle=True)
heat3 = np.load('heat3.npy', allow_pickle=True)
heat_tot = ((heat/len1)**2 + (heat2/len2)**2 + (heat3/len3)**2)**(1/2)


def pca(heat) :
    heat_std = StandardScaler().fit_transform(heat)
    sklearn_pca = sklearnPCA(n_components = len(heat))
    Y_sklearn = np.round(sklearn_pca.fit_transform(heat_std), 3)
    percentage = np.round(sklearn_pca.explained_variance_[0:2]/np.sum(sklearn_pca.explained_variance_),3)
    return Y_sklearn[:,0:2],percentage

def plot_pca(Y_sklearn, name, percentage) :
    with plt.style.context('seaborn-whitegrid'):
        plt.figure(figsize=(8, 6))
        plt.scatter(Y_sklearn[:,0],Y_sklearn[:,1])
        plt.xlabel('Principal Component 1 : ' + str(round(percentage[0]*100,3)) + '%')
        plt.ylabel('Principal Component 2 : ' + str(round(percentage[1]*100,3)) + '%')
        plt.legend(loc='lower center')
        plt.tight_layout()
        plt.savefig(name + '.png', format='png', dpi=1200)
        plt.show()
        
        
def kmeans_cluster(Y_sklearn, n, name, percentage) :
    kmeans = KMeans(n_clusters=n, random_state=0).fit(Y_sklearn)
    clusters = kmeans.labels_
    df = pd.DataFrame(dict(x=Y_sklearn[:,0], y=Y_sklearn[:,1], label=clusters))
    colors = {0:'red', 1:'blue', 2:'green', 3:'orange', 4:'gold', 5:'purple', 6:'grey', 7:'pink',\
              8:'navy', 9:'springgreen', 10:'salmon', 11:'skyblue', 12:'tan', 13:'sienna',\
              14:'turquoise', 15:'teal', 16:'chartreuse', 17:'crimson', 18:'fuchsia', 19:'beige',\
              20:'yellow', 21:'aqua', 22:'olivedrab', 23:'deeppink', 24:'maroon', 25:'mistyrose',\
              26:'seagreen', 27:'darkorange', 28:'mediumpurple', 29:'khaki'}
    plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots()
    grouped = df.groupby('label')
    for key, group in grouped:
        group.plot(ax=ax, kind='scatter', x='x', y='y', label=key, color=colors[key])
    plt.legend(prop={'size':8})
    plt.rc('xtick', labelsize=6) 
    plt.rc('ytick', labelsize=6) 
    plt.xlabel('Principal Component 1 : ' + str(round(percentage[0]*100,3)) + '%', fontsize = 10)
    plt.ylabel('Principal Component 2 : ' + str(round(percentage[1]*100,3)) + '%', fontsize = 10)
    plt.savefig(name + '.png', format ='png', dpi=1200)
    plt.close()
    return clusters


def agglomerative_cluster(Y_sklearn, n, norm, link, name, percentage) :
    agglomerative = AgglomerativeClustering(n_clusters = n, affinity = norm, compute_full_tree = 'auto', linkage = link).fit(Y_sklearn)
    clusters = agglomerative.labels_
    df = pd.DataFrame(dict(x=Y_sklearn[:,0], y=Y_sklearn[:,1], label=clusters))
    colors = {0:'red', 1:'blue', 2:'green', 3:'orange', 4:'gold', 5:'purple', 6:'grey', 7:'pink',\
              8:'navy', 9:'springgreen', 10:'salmon', 11:'skyblue', 12:'tan', 13:'sienna',\
              14:'turquoise', 15:'teal', 16:'chartreuse', 17:'crimson', 18:'fuchsia', 19:'beige',\
              20:'yellow', 21:'aqua', 22:'olivedrab', 23:'deeppink', 24:'maroon', 25:'mistyrose',\
              26:'seagreen', 27:'darkorange', 28:'mediumpurple', 29:'khaki'}
    plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots()
    grouped = df.groupby('label')
    for key, group in grouped:
        group.plot(ax=ax, kind='scatter', x='x', y='y', label=key, color=colors[key])
    plt.legend(prop={'size':8})
    plt.rc('xtick', labelsize=6) 
    plt.rc('ytick', labelsize=6) 
    plt.xlabel('Principal Component 1 : ' + str(round(percentage[0]*100,3)) + '%', fontsize = 10)
    plt.ylabel('Principal Component 2 : ' + str(round(percentage[1]*100,3)) + '%', fontsize = 10)
    plt.savefig(name + '.png', format ='png', dpi=1200)
    plt.close()
    return clusters



def top_all(multi, n, n_top) :
    multi_top = [list() for i in range(len(n_top))]
    multi_topL = [list() for i in range(len(n_top))]
    multi_topL_val = list()
    for k in range(len(n_top)) :
        for i in range(len(multi)) : 
            multi_top[k].append(multi[i]['node'][0:n_top[k]])
        for i in range(len(multi)) : 
            multi_topL[k].append(list(multi_top[k][i])) 
    for k in range(len(n_top)) :
        multi_topL[k] = list(itertools.chain(*multi_topL[k]))
    for k in range(len(n_top)) :
        multi_topL_val.append( Counter(multi_topL[k]))
    return multi_top, multi_topL_val
    


def top_cluster(multi_top, cluster, n, n_top) :
    multi_topL = [list() for i in range(n)]
    multi_topL_val = [list() for i in range(n)]
    index_cluster = list()
    multi_sub = [list() for i in range(n)]
    for k in range(n) :
        index_cluster.append(np.where(cluster == k))
        multi_sub[k] = [list() for i in range(len(n_top))]
        multi_topL[k] = [list() for i in range(len(n_top))]
        multi_topL_val[k] = [list() for i in range(len(n_top))]
    for k in range(n) :
        for l in range(len(n_top)) :
            multi_sub[k][l].append([multi_top[l][i] for i in index_cluster[k][0]])
    for k in range(n) :
        for l in range(len(n_top)) :
            for i in range(len(index_cluster[k][0])) :
                multi_topL[k][l].append(list(multi_sub[k][l][0][i]))
    for k in range(n) :
        for l in range(len(n_top)) :
            multi_topL[k][l] = list(itertools.chain(*multi_topL[k][l]))
    for k in range(n) :
        for l in range(len(n_top)) :         
            multi_topL_val[k][l].append(Counter(multi_topL[k][l]))
    return multi_sub, multi_topL_val


    
def plot_bar_cluster(multi_cluster_topL, n, n_tot, name, cluster) :
    count_cluster = Counter(cluster)
    for k in range(len(n_tot)) :
        size = int(n/2) + n%2
        fig,a =  plt.subplots(size,2)
        plt.subplots_adjust(hspace = 0.75)
        compt = 0
        compt2 = 0
        for l in range(n) :
            keys = list(multi_cluster_topL[l][k][0].keys())
            values = list(multi_cluster_topL[l][k][0].values())
            y_pos = np.arange(len(keys))
            barlist = a[compt][compt2].bar(y_pos, values, tick_label = keys)
            max_val = count_cluster[l]
            index_max = np.where(np.array(values) == max_val)[0]
            for i in range(len(index_max)) :
                barlist[index_max[i]].set_color('r')
            a[compt][compt2].set_title('cluster ' + str(l), size = 6)
            a[compt][compt2].set_xticklabels([])
            a[compt][compt2].tick_params(axis = 'x', which = 'both', bottom = False,\
                    top  = False, labelbottom = False)
            compt2 += 1
            if (compt2%2 == 0):
                compt2 = 0
                compt += 1
        plt.savefig(name + '_Distribution_Top' + str(n_tot[k]) + '.png', dpi=1200)
        
def plot_bar(multi_topL, n_tot, name) :
    size = int(len(n_tot)/2) + len(n_tot)%2
    fig,a =  plt.subplots(size,2)
    plt.subplots_adjust(hspace = 0.5)
    compt = 0
    compt2 = 0
    for k in range(len(n_tot)) :
        keys = list(multi_topL[k].keys())
        values = list(multi_topL[k].values())
        y_pos = np.arange(len(keys))
        barlist = a[compt][compt2].bar(y_pos, values, tick_label = keys)
        max_val = num
        index_max = np.where(np.array(values) == max_val)[0]
        for i in range(len(index_max)) :
            barlist[index_max[i]].set_color('r')
        a[compt][compt2].set_title('Top' + str(n_tot[k]))
        a[compt][compt2].set_xticklabels([])
        a[compt][compt2].tick_params(axis = 'x', which = 'both', bottom = False,\
                    top  = False, labelbottom = False)
        compt2 += 1
        if (compt2%2 == 0):
            compt2 = 0
            compt += 1
    plt.savefig(name + '_Distribution_all.png', dpi=1200)


n_top = [1, 3, 5, 10, 20, 100]
n_clus1 = 4
Y_sklearn1, per1 = pca(heat)
plot_pca(Y_sklearn1, 'heat_pca', per1)
clusters1 = agglomerative_cluster(Y_sklearn1, n_clus1, 'l1', 'average', 'heat_agglomerative', per1)
clusters1 = kmeans_cluster(Y_sklearn1, n_clus1, 'heat_kmeans', per1)
multi1_top, multi1_topL_val = top_all(multi_1, n_clus1, n_top)
multi1_sub, multi1_cluster_topL_val = top_cluster(multi1_top, clusters1, n_clus1, n_top)
plot_bar(multi1_topL_val, n_top, 'heat')
plot_bar_cluster(multi1_cluster_topL_val, n_clus1, n_top,'heat', clusters1)

n_clus2 = 7
Y_sklearn2, per2 = pca(heat2)
plot_pca(Y_sklearn2, 'heat2_pca', per2)
clusters2 = agglomerative_cluster(Y_sklearn2, n_clus2, 'l1', 'average', 'heat2_agglomerative', per2)
clusters2 = kmeans_cluster(Y_sklearn2, n_clus2, 'heat2_kmeans', per2)
multi2_top, multi2_topL_val = top_all(multi_2, n_clus2, n_top)
multi2_sub, multi2_cluster_topL_val = top_cluster(multi2_top, clusters2, n_clus2, n_top)
plot_bar(multi2_topL_val, n_top, 'heat2')
plot_bar_cluster(multi2_cluster_topL_val, n_clus2, n_top,'heat2', clusters2)

n_clus3 = 4
Y_sklearn3, per3 = pca(heat3)
plot_pca(Y_sklearn3, 'heat3_pca', per3)
clusters3 = agglomerative_cluster(Y_sklearn3, n_clus3, 'l1', 'average', 'heat3_agglomerative', per3)
clusters3 = kmeans_cluster(Y_sklearn3, n_clus3, 'heat3_kmeans', per3)
multi3_top, multi3_topL_val = top_all(multi_3, n_clus3, n_top)
multi3_sub, multi3_cluster_topL_val = top_cluster(multi3_top, clusters3, n_clus3, n_top)
plot_bar(multi3_topL_val, n_top, 'heat3')
plot_bar_cluster(multi3_cluster_topL_val, n_clus3, n_top,'heat3', clusters3)



n_clus_tot = 6
Y_sklearn_tot, per_tot = pca(heat_tot)
plot_pca(Y_sklearn_tot, 'heattot_pca', per_tot)
clusters_tot = agglomerative_cluster(Y_sklearn_tot, n_clus_tot, 'l1', 'average', 'heattot_agglomerative', per_tot)
clusters_tot = kmeans_cluster(Y_sklearn_tot, n_clus_tot, 'heattot_kmeans', per_tot)

multi1tot_top, multi1tot_topL_val = top_all(multi_1, n_clus_tot, n_top)
multi1tot_sub, multi1tot_cluster_topL_val = top_cluster(multi1_top, clusters_tot, n_clus_tot, n_top)
plot_bar(multi1tot_topL_val, n_top, 'heattot')
plot_bar_cluster(multi1tot_cluster_topL_val, n_clus_tot, n_top,'heattot', clusters_tot)

multi2tot_top, multi2tot_topL_val = top_all(multi_2, n_clus_tot, n_top)
multi2tot_sub, multi2tot_cluster_topL_val = top_cluster(multi2tot_top, clusters_tot, n_clus_tot, n_top)
plot_bar(multi2tot_topL_val, n_top, 'heat2tot')
plot_bar_cluster(multi2tot_cluster_topL_val, n_clus_tot, n_top,'heat2tot', clusters_tot)

multi3tot_top, multi3tot_topL_val = top_all(multi_3, n_clus_tot, n_top)
multi3tot_sub, multi3tot_cluster_topL_val = top_cluster(multi3tot_top, clusters_tot, n_clus_tot, n_top)
plot_bar(multi3tot_topL_val, n_top, 'heat3tot')
plot_bar_cluster(multi3tot_cluster_topL_val, n_clus_tot, n_top,'heat3tot', clusters_tot)

    
