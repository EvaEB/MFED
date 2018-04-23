#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:42:06 2017

will generate new parameters based on distances found in all stats files
(located in [stats_dir] and names stats_{something}.npy) and [max_dist].
New parameter files are stored in [new_par_dir] and named par_set_[name].npy

usage: python next_par.py [max_dist] [stats_dir] [new_par_dir] [name]
"""
import numpy as np
from collections import Counter
import scipy.stats as stat
import sys
import os

def in_prior(par,model):
    if model == 0: #beneficials only
        if par[0]<0 or par[0]>1:
            return False
        if par[1]<0 or par[1]>2:
            return False
        if par[2]<-10 or par[2]>np.log(160):
            return False
        return True
    if model == 1: #lethals only
        if par[0]<0 or par[0]>1:
            return False
        if par[1]<-10 or par[1]>np.log(160):
            return Falsene
        return True
    if model == 2: #lethals and beneficials
        if par[0]<0 or par[0]>1:
            return False
        if par[1]<0 or par[1]>1:
            return False
        if par[2]<0 or par[2]>2:
            return False
        if par[0]+par[1]>1:
            return False
        if par[3]<-10 or par[3]>np.log(160):
            return False
        else:
            return True
    if model == 3 or model == 4: #lognormal (truncated)
        if par[0]<0 or par[0]>1:
            return False
        if par[1]<-1 or par[1]>1:
            return False
        if par[2]<=0 or par[2]>1:
            return False
        if par[3]<-10 or par[3]>np.log(160):
            return False
        else:
            return True
    if model == 5: #spikes
        if par[0]<0 or par[0]>1:
            return False
        if par[1]<0 or par[1]>1:
            return False
        if par[2]<0 or par[2]>1:
            return False
        if par[3]<0 or par[3]>1:
            return False
        if par[4]<0 or par[4]>1:
            return False
        if par[5]<0 or par[5]>1:
            return False
        if par[6]<1 or par[6]>2:
            return False
        if par[7]<-10 or par[7]>np.log(160):
            return False
        if par[0]+par[1]+par[2]+par[3]>1:
            return False
        if par[4]>par[5]:
            return False
        else:
            return True
    if model == 6: #neutral
        if par[0]<-10 or par[0]>np.log(160):
            return False
        return True

def get_weights(par,old_par,s,weights,model_prop):
    kernel = [np.prod([stat.norm.pdf(par[i],old_par[j,i],s[i]) for i in range(len(par))]) for j in range(len(old_par))]
    return 1/(sum(kernel*weights)*model_prop)




max_sims = 100
max_dist = float(sys.argv[1])
n_per_iter = 500
simdir = sys.argv[2]
new_dir = sys.argv[3]
name = sys.argv[4]

n_par = [3,2,4,4,4,8,1]
first=True
new_par_set = []

for i in range(1,max_sims):
    if first:
        try:
            stats = np.load('{}/stats_{}.npy'.format(simdir,i))
            first=False
        except IOError:
            first=True
            print 'not found: {}'.format(i)
    else:
        try:
            stats = np.vstack((stats, np.load('{}/stats_{}.npy'.format(simdir,i))))
        except IOError:
            print 'not found: {}'.format(i)
accepted = stats[stats[:,-1]<max_dist]
models = np.random.choice(accepted[:,0],n_per_iter)
models = Counter(models)
for i in models:
    #get the model
    current_model = accepted[accepted[:,0]==i]

    #get the sigmas
    sigmas = []
    i = int(i)
    for j in range(n_par[i]):
        sigmas.append(np.std(current_model[:,j+1])/5)
        if sigmas[-1] == 0:
            sigmas[-1] = np.std(stats[(stats[:,0]==i)][:,j+1])/5

    #normalize the weights
    weights = current_model[:,9]/sum(current_model[:,9])
    #choose the parameters to continue with
    chosen = np.random.choice(range(len(current_model)),size=models[i],p=weights)
    for par_set in chosen:
        pars = current_model[par_set,1:n_par[i]+1]
        if pars[-1] > 0:
            pars[-1] = np.log(pars[-1])
        while True:
            new_pars = []
            for j, par in enumerate(pars):
                new_pars.append(np.random.normal(par,sigmas[j]))
            if in_prior(new_pars,i):
                break
        weight = get_weights(new_pars,current_model[:,1:n_par[i]+1],sigmas,weights,1)
        padding = [0 for k in range(8-len(pars))]


        new_par_set.append([i]+new_pars+padding+[weight])
new_par_set = np.array(new_par_set)
np.save('{}/par_{}.npy'.format(new_dir,name),new_par_set)
