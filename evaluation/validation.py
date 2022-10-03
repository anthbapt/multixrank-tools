#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import multixrank as mxk
import multiprocessing as mp
import pandas as pd
import copy
import glob
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)

try : 
    os.mkdir('seeds')
    os.mkdir('parameters')
    os.mkdir('bipartites')
    os.mkdir('results')
except OSError : 
    pass
    

def mxrank(k) :
    print(str(k) + '\n')
    multixrank_obj = mxk.Multixrank(config = "parameters/config_full_" + str(k) + ".yml", wdir = path)
    ranking = multixrank_obj.random_walk_rank()
    multixrank_obj.write_ranking(ranking, path = "results/seeds_" + str(k))

def create_parameters(k, path) : 
    r = 0.7
    eta = [0.5, 0.5, 0]
    lamb = np.array([[1/3,1/3,1/3],
            [1/3,1/3,1/3],
            [1/3,1/3,1/3]])
    delta1 = 0.5
    delta2 = 0
    delta3 = 0.5
    tau = [[1/3, 1/3, 1/3],
           [1],
           [1/4, 1/4, 1/4, 1/4]]
    file = open(path + 'parameters/' + 'config_full_' + str(k) + '.yml','w')
    file.write('seed: seeds/seeds_' + str(k) + '.txt' + '\n')
    file.write('self_loops: 0' + '\n')
    file.write('r: ' + str(r) + '\n')
    file.write('eta: ' + '[{},{},{}]'.format(eta[0],eta[1],eta[2]) + '\n')
    file.write('lamb:' + '\n' + '    ' + \
               '- [{},{},{}]'.format(lamb[0,0], lamb[0,1], lamb[0,2]) + '\n' + '    '  + \
               '- [{},{},{}]'.format(lamb[1,0], lamb[1,1], lamb[1,2]) + '\n' + '    '  + \
               '- [{},{},{}]'.format(lamb[2,0], lamb[2,1], lamb[2,2]) + '\n')
    file.write('multiplex:' + '\n' + '    ')
    file.write('protein:' + '\n' + '        ' + \
                       'layers:' + '\n' + '            ' + \
                           '- networks/multiplex/1/PPI.gr' + '\n' + '            ' + \
                           '- networks/multiplex/1/Complexes.gr' + '\n' + '            ' + \
                           '- networks/multiplex/1/Reactome.gr' + '\n' + '        ' + \
                        'delta: {}'.format(delta1) + '\n' + '        ' + \
                        'graph_type: ' + '[00, 00, 00]' + '\n' + '        ' + \
                        'tau: ' + str(tau[0]) + '\n')
    file.write('    ' + 'disease:' + '\n' + '        ' + \
                       'layers:' + '\n' + '            ' + \
                           '- networks/multiplex/2/disease_disease_final.gr' + '\n' + '        ' + \
                        'delta: {}'.format(delta2) + '\n' + '        ' + \
                        'graph_type: ' + '[00]' + '\n' + '        ' + \
                        'tau: ' + str(tau[1]) + '\n')
        
    file.write('    ' + 'drug:' + '\n' + '        ' + \
                       'layers:' + '\n' + '            ' + \
                           '- networks/multiplex/3/clinical_drug_interactions.gr' + '\n' + '            ' + \
                           '- networks/multiplex/3/experimental_drug_combinations.gr' + '\n' + '            ' + \
                           '- networks/multiplex/3/predicted_drug_combinations.gr' + '\n' + '            ' + \
                           '- networks/multiplex/3/drugbank-chem-chem.gr' + '\n' + '        ' + \
                        'delta: {}'.format(delta3) + '\n' + '        ' + \
                        'graph_type: ' + '[00, 00, 00, 00]' + '\n' + '        ' + \
                        'tau: ' + str(tau[2]) + '\n')
    file.write('bipartite: ' + '\n' + '    ' + \
               'bipartites/1_2_' + str(k) + '.gr: ' + '{source: protein, target: disease, graph_type: 00}' + '\n' + '    ' + \
               'networks/bipartite/3_1.gr: ' + '{source: drug, target: protein, graph_type: 00}' + '\n' + '    ' + \
               'networks/bipartite/3_2.gr: ' + '{source: drug, target: disease, graph_type: 00}')
    file.close
    
    
def protocol1(path, nodes_1, bipartite, bip, num_cpu) :
    os.chdir(path)
    results = pd.DataFrame(columns = ['node1', 'node2', 'rank'])
    compt = 0
    list_nodes = list()
    for i in nodes_1 :
        if (bip == '1_2') :
            associated_node = set(bipartite[bipartite[1] == i][0])
        else :
            associated_node = set(bipartite[bipartite[0] == i][1])
        for j in associated_node :
            temp = pd.read_csv(path + bip + '.gr', sep ='\t', header = None)
            if (bip == '1_2') :
                index = temp[(temp[0] == j) & (temp[1] == i)].index
                new_bipartite = temp.drop(index)
                reverse_associated_node = set(bipartite[bipartite[0] == j][1])
                for k in reverse_associated_node :
                    index = new_bipartite[(new_bipartite[0] == j) & (new_bipartite[1] == k)].index
                    new_bipartite2 = new_bipartite.drop(index)
                list_nodes.append((i,j))
            else : # bip == '2_1'
                index = temp[(temp[1] == j) & (temp[0] == i)].index
                new_bipartite = temp.drop(index)
                reverse_associated_node = set(bipartite[bipartite[1] == j][0])
                for k in reverse_associated_node :
                    index = new_bipartite[(new_bipartite[1] == j) & (new_bipartite[0] == k)].index
                    new_bipartite2 = new_bipartite.drop(index) 
                list_nodes.append((j,i))
            new_bipartite2.to_csv('bipartites/' + bip + '_' + str(compt) + '.gr', sep = '\t', header = None, index = False)
            seed = str(i)
            seed_file = open('seeds/seeds_' + str(compt) + '.txt', 'w')
            seed_file.write(seed)
            seed_file.close()
            compt += 1
    size = len(glob.glob("seeds/*.txt"))
    for k in range(size) :
        create_parameters(k, path)
    p = mp.Pool(processes=num_cpu)
    p.map(mxrank, [i for i in range(size)])    
    for k in range(size) : 
        ranking = pd.read_csv('results/seeds_' + str(k) + '/multiplex_protein.tsv', sep = '\t')
        results = results.append({'node1' : list_nodes[k][0], 'node2' : list_nodes[k][1], 'rank' : ranking[ranking['node'] == list_nodes[k][1]].index[0]}, ignore_index=True)
    return results
    

def protocol2(path, nodes_1, bipartite, bip, num_cpu) :
    os.chdir(path)
    results = pd.DataFrame(columns = ['node1', 'node2', 'rank'])
    compt = 0
    list_nodes = list()
    for i in nodes_1 :
        if (bip == '1_2') :
            associated_node = set(bipartite[bipartite[1] == i][0])
        else :
            associated_node = set(bipartite[bipartite[0] == i][1])
        for j in associated_node :
            temp = pd.read_csv(path + bip + '.gr', sep ='\t', header = None)
            if (bip == '1_2') :
                index = temp[(temp[0] == j) & (temp[1] == i)].index
                new_bipartite = temp.drop(index)
                list_nodes.append((i,j))
            else : # bip == '2_1'
                index = temp[(temp[1] == j) & (temp[0] == i)].index
                new_bipartite = temp.drop(index)
                list_nodes.append((j,i))
            new_bipartite.to_csv('bipartites/' + bip + '_' + str(compt) + '.gr', sep = '\t', header = None, index = False)
            seed = str(i)
            seed_file = open('seeds/seeds_' + str(compt) + '.txt', 'w')
            seed_file.write(seed)
            seed_file.close()
            compt += 1
    size = len(glob.glob("seeds/*.txt"))
    for k in range(size) :
        create_parameters(k, path)
    p = mp.Pool(processes=num_cpu)
    p.map(mxrank, [i for i in range(size)])    
    for k in range(size) : 
        ranking = pd.read_csv('results/seeds_' + str(k) + '/multiplex_protein.tsv', sep = '\t')
        results = results.append({'node1' : list_nodes[k][0], 'node2' : list_nodes[k][1], 'rank' : ranking[ranking['node'] == list_nodes[k][1]].index[0]}, ignore_index=True)
    return results


def loocv(path, nodes_1, bipartite, bip, num_cpu):
    os.chdir(path)
    results = pd.DataFrame(columns = ['node1', 'node2', 'rank'])
    compt = 0
    list_nodes = list()
    for i in nodes_1 :
        if (bip == '1_2') :
            associated_node = set(bipartite[bipartite[1] == i][0])
        else :
            associated_node = set(bipartite[bipartite[0] == i][1])
        if (len(associated_node)) >= 2 :
            for j in associated_node :
                temp = pd.read_csv(path + bip + '.gr', sep ='\t', header = None)
                if (bip == '1_2') :
                    index = temp[(temp[0] == j) & (temp[1] == i)].index
                    new_bipartite = temp.drop(index)
                    list_nodes.append((i,j))
                else : # bip == '2_1'
                    index = temp[(temp[1] == j) & (temp[0] == i)].index
                    new_bipartite = temp.drop(index)
                    list_nodes.append((j,i))
                new_bipartite.to_csv('bipartites/' + bip + '_' + str(compt) + '.gr', sep = '\t', header = None, index = False)
                seeds = list()
                seeds.append(i)
                temp = copy.deepcopy(associated_node)
                temp.remove(j)
                for k in temp :
                    seeds.append(k)
                seed_file = open('seeds/seeds_' + str(compt) + '.txt', 'w')
                for k in seeds :
                    seed_file.write(str(k) + '\n')
                seed_file.close()
                compt += 1
    size = len(glob.glob("seeds/*.txt"))
    for k in range(size) :
        create_parameters(k, path)
    p = mp.Pool(processes=num_cpu)
    p.map(mxrank, [i for i in range(size)])    
    for k in range(size) : 
        ranking = pd.read_csv('results/seeds_' + str(k) + '/multiplex_protein.tsv', sep = '\t')
        results = results.append({'node1' : list_nodes[k][0], 'node2' : list_nodes[k][1], 'rank' : ranking[ranking['node'] == list_nodes[k][1]].index[0]}, ignore_index=True)
    return results
