import multixrank as mxk
import multiprocessing as mp
import glob
import os
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)

try : 
    os.mkdir('results')
except OSError : 
    pass
    

def mxrank(k) :
    print(str(k) + '\n')
    multixrank_obj = mxk.Multixrank(config = "parameters/config_full_" + str(k) + ".yml", wdir = path)
    ranking = multixrank_obj.random_walk_rank()
    multixrank_obj.write_ranking(ranking, path = "results/param_" + str(k))


size = len(glob.glob("parameters/*.yml"))
num_cpu = 14
p = mp.Pool(processes=num_cpu)
p.map(mxrank, [i for i in range(size)])    
    

