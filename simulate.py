# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 17:28:24 2016

simulation of sequence evolution

either import 'simulate' into a script, or run simulate.py -h
for info on arguments 
"""
import numpy as np
import random
import sys
import os
from simulation_functions import *
import argparse


def simulate(patients, simdir, model, par,save_fitness=False):
    N = 2600

    f = open(patients)
    settings = f.readlines()
#    f.close()
#    with open(simdir+ '\\test.txt', 'w') as f:
#        f.write(str(settings))

    nrMut = np.array(settings[1].strip().split(),dtype='int')
    samples = np.array(settings[3].strip().split(),dtype='int')
    GAincrease = np.array(settings[5].strip().split(),dtype='int')
    ind = len(samples)

    seq = createSeq(N)
    fitnessTable = createFitnessTable(seq, model, par[:-1])
    if save_fitness:
        np.save(simdir+'/fitnessTable.npy',fitnessTable)
    for ID in range(len(samples)):
        if GAincrease[ID] == 1:
            simulate_and_print(seq,fitnessTable,nrMut[ID],samples[ID],ID,simdir,GAincrease=par[-1])
        else:
            simulate_and_print(seq,fitnessTable,nrMut[ID],samples[ID],ID,simdir,GAincrease=1)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='simulate a dataset')
    parser.add_argument('-patients',help='path to the patients file')
    parser.add_argument('-simdir',help='path to directory for simulation output')
    parser.add_argument('-model',help='model for the MFED. one of "neutral",\
                                       "lethals only","beneficials only",\
                                       "lethals and beneficials","lognormal",\
                                       "lognormal truncated","spikes"')
    parser.add_argument('-par',nargs=argparse.REMAINDER,
                        help='the parameters of the MFED')
    parser.add_argument('-save',action='store_true', default=False,
                        help='save the fitnessTable as a .npy file')


    args = parser.parse_args()

    if args.model == 'lethals only':
        try:
            par = [args.par[args.par.index('fl')+1],
                   args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'fl or GA not provided in par, exiting'
            exit()
    elif args.model == 'beneficials only':
        try:
            par = [args.par[args.par.index('fb')+1],
                   args.par[args.par.index('beta')+1],
                   args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'fb, beta or GA not provided in par, exiting'
            exit()
    elif args.model == 'lethals and beneficials':
        try:
            par = [args.par[args.par.index('fl')+1],
                   args.par[args.par.index('fb')+1],
                   args.par[args.par.index('beta')+1],
                   args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'fl, fb, beta or GA not provided in par, exiting'
            exit()
    elif 'lognormal' in args.model:
        try:
            par = [args.par[args.par.index('fl')+1],
                   args.par[args.par.index('mu')+1],
                   args.par[args.par.index('sigma')+1],
                   args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'fl, mu, sigma or GA not provided in par, exiting'
            exit()
    elif args.model == 'neutral':
        try:
            par = [args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'GA not provided in par, exiting'
            exit()
    elif '2exp' in args.model:
        try:
            par = [args.par[args.par.index('fd')+1],
                   args.par[args.par.index('fb')+1],
                   args.par[args.par.index('ld')+1],
                   args.par[args.par.index('lb')+1],
                   args.par[args.par.index('GA')+1]]
        except ValueError:
            print 'f_(d/b), l(d/b) or GA not provided in par, exiting'
            exit()
    else:
        print 'argument parsing for {} not yet implemented or unknown model, \
               exiting'.format(args.model)
        exit()

    par = [float(i) for i in par]
    par[-1] = np.exp(par[-1]) #transform G-A rate (provided as log)
    print par[-1]
    simulate(args.patients, args.simdir,args.model,par,save_fitness=args.save)
