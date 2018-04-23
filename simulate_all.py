#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:41:26 2017

performs a full simulation for all parameter sets in [parameter_file] and
calculates their summary statistics. Temporary files are stored in
[simulation_directory].
output is reported in [outfile], a .npy file (binary file of numpy array), with
a row for each given parameter set, and 25 columns
1st column: model ID
2-9th column: parameters, padded with zeros
10th column: weight of parameter sets
11-24th column: summary statistics of simulated dataset
25th column: distance to the data

usage: python simulate_all.py [parameter_file] [simulation_directory] [outfile]

"""

import numpy as np
import simulate
import sys
import os
import getSummaryStats

def distance(data_sumstat, sumstats):
    return sum((np.ones(len(data_sumstat))-(np.array(sumstats)/np.array(data_sumstat)))**2)

par_file = sys.argv[1]
patients_file = 'patients'
simdir = sys.argv[2]
outfile = sys.argv[3]

sumstats_data = [0.830543754873,0.741666666667,0.0888770882064,0.311894779355,
                 963,164,71,8,1,0.46945061741,0.531440032429,0.0148259877411,
                 0.635384983705,0.404868355688]


par = np.load(par_file)

models = ['beneficials only','lethals only','lethals and beneficials',
          'lognormal','lognormal truncated','spikes','neutral']

n_par  = [3,                    2,              4,
          4,                   4,               8,      1]

#check progress on outfile
try:
    all_stats = list(np.load(outfile))
    start = len(all_stats)
    print 'restarting from',start
except IOError:
    start = 0
    all_stats = []


for i in par[start:]:
    model = models[int(i[0])]
    par = i[1:n_par[int(i[0])]+1]
    par[-1] = np.exp(par[-1])
    print par[-1]
    simulate.simulate(patients_file,simdir,model, par)

    #move all simulations to one file
    filenames = os.listdir(simdir)
    with open(simdir+'/out', 'w') as fout:
        for name in filenames:
            if name != 'out':
                with open(simdir+'/'+name) as fin:
                    for line in fin.readlines():
                        fout.write(line)
                os.remove(simdir+'/'+name)

    #get summary statistics
    stats=getSummaryStats.get_stats(simdir+'/out')
    dist = distance(sumstats_data,stats)
    all_stats.append(list(i)+stats+[dist])
    np.save(outfile,np.array(all_stats))
