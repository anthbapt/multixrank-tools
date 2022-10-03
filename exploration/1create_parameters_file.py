#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import os
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)

try : 
    os.mkdir('parameters')
except OSError : 
    pass


r = 0.7
tau = [[1/3,1/3,1/3],\
       [1],\
       [1/4,1/4,1/4,1/4]]
delta1 = np.linspace(0.1,1,3)
delta2 = np.linspace(0.1,1,3)
delta3 = np.linspace(0.1,1,3)
lambda21 = np.linspace(0.1,0.9,3)
lambda31 = np.linspace(0.1,0.9,3)
lambda11 = 1-lambda21-lambda31
lambda12 = lambda21
lambda32 = np.linspace(0.1,0.9,3)
lambda22 = 1-lambda12-lambda32
lambda13 = lambda31
lambda23 = lambda32
lambda33 = 1-lambda13-lambda23
eta1 = 1
eta2 = 0
eta3 = 0
compt = 0

for b in delta1 :
    for c in delta2 :
        for a in delta3 :
            for d in lambda12 :
                for e in lambda13 :
                    for f in lambda23 :
                        if (d+e<1) and (d+f)<1 and (e+f)<1 :
                            print(compt,b,c,a,d,e,f)
                            file = open(path + 'parameters/' + 'config_full_' + str(compt) + '.yml','w')
                            file.write('seed: seeds.txt' + '\n')
                            file.write('self_loops: 0' + '\n')
                            file.write('r: ' + str(r) + '\n')
                            file.write('eta: ' + '[{},{},{}]'.format(eta1,eta2,eta3) + '\n')
                            file.write('lamb:' + '\n' + '    ' + \
                                       '- [{},{},{}]'.format(1-d-e, d, e) + '\n' + '    '  + \
                                       '- [{},{},{}]'.format(d, 1-d-f, f) + '\n' + '    '  + \
                                       '- [{},{},{}]'.format(e, f, 1-e-f) + '\n')
                            file.write('multiplex:' + '\n' + '    ')
                            file.write('protein:' + '\n' + '        ' + \
                                               'layers:' + '\n' + '            ' + \
                                                   '- networks/multiplex/1/PPI.gr' + '\n' + '            ' + \
                                                   '- networks/multiplex/1/Complexes.gr' + '\n' + '            ' + \
                                                   '- networks/multiplex/1/Reactome.gr' + '\n' + '        ' + \
                                                'delta: {}'.format(b) + '\n' + '        ' + \
                                                'graph_type: ' + '[00, 00, 00]' + '\n' + '        ' + \
                                                'tau: ' + str(tau[0]) + '\n')
                            file.write('    ' + 'disease:' + '\n' + '        ' + \
                                               'layers:' + '\n' + '            ' + \
                                                   '- networks/multiplex/2/disease_disease_final.gr' + '\n' + '        ' + \
                                                'delta: {}'.format(c) + '\n' + '        ' + \
                                                'graph_type: ' + '[00]' + '\n' + '        ' + \
                                                'tau: ' + str(tau[1]) + '\n')
                                
                            file.write('    ' + 'drug:' + '\n' + '        ' + \
                                               'layers:' + '\n' + '            ' + \
                                                   '- networks/multiplex/3/clinical_drug_interactions.gr' + '\n' + '            ' + \
                                                   '- networks/multiplex/3/experimental_drug_combinations.gr' + '\n' + '            ' + \
                                                   '- networks/multiplex/3/predicted_drug_combinations.gr' + '\n' + '            ' + \
                                                   '- networks/multiplex/3/drugbank-chem-chem.gr' + '\n' + '        ' + \
                                                'delta: {}'.format(a) + '\n' + '        ' + \
                                                'graph_type: ' + '[00, 00, 00, 00]' + '\n' + '        ' + \
                                                'tau: ' + str(tau[2]) + '\n')
                            file.write('bipartite: ' + '\n' + '    ' + \
                                       'networks/bipartite/2_1.gr: ' + '{source: disease, target: protein, graph_type: 01}' + '\n' + '    ' + \
                                       'networks/bipartite/3_1.gr: ' + '{source: drug, target: protein, graph_type: 00}' + '\n' + '    ' + \
                                       'networks/bipartite/3_2.gr: ' + '{source: drug, target: disease, graph_type: 00}')
                            file.close
                            compt += 1
