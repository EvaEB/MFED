#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:12:49 2017

generates 100 random parameters for all 7 defined models
(0: beneficials only, 1: lethals only, 2: lethals and beneficials,
3: lognormal, 4: lognormal truncated, 5: 5 spikes and 6: neutral)

usage: python initial_smc.py [filename]

output: a .npy file (binary file of numpy array), with 700 rows (100 per model)
and 10 rows
model ID - parameters (padded with 0s if there's less than 8) - weight (= 1)
"""
import numpy as np
import sys

def simulation_par(model,prior,N,maxPar,prior_condition):
    filler = [0 for i in range(maxPar-len(prior))]
    all_pars = []
    for i in range(N):
        while True:
            pars = [np.random.uniform(prior[i][0],prior[i][1]) for i in range(len(prior))]
            dist = [par[2] for par in prior]
            if 'log' in dist:
                logs = [k for k in range(len(dist)) if dist[k] == 'log']
                for par in logs:
                    pars[par] = np.exp(pars[par])
            if prior_condition(pars):
                break
        all_pars.append([model]+pars+filler+[1])
    return np.array(all_pars)

#model definitions
#0  beneficials only
#1  lethals only
#2  beneficials and lethals
#3  log-normal
#4  log-normal truncated
#5  5-spikes
#6  neutral

N = 100   #samples per model
prior = {}
maxPar = 8

models = [0,1,2,3,4,5,6]
ga_prior = [-10,np.log(160),'uniform']

#0  beneficials only
if 0 in models:
    def prior_condition_none(pars):
        return True

    model = 0
    prior = [[0,1,'uniform'],[0,2,'uniform'],ga_prior]

    par = simulation_par(model,prior,N,maxPar,prior_condition_none)

#1  lethals only
if 1 in models:
    model = 1
    prior = [[0,1,'uniform'],ga_prior]

    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_none)))

#2  beneficials and lethals
if 2 in models:
    def prior_condition_ben_lethal(pars):
        if pars[0]+pars[1]<=1:
            return True
        return False


    model = 2
    prior = [[0,1,'uniform'],[0,1,'uniform'],[0,2,'uniform'],ga_prior]

    #par = simulation_par(model,prior,N,maxPar,prior_condition_ben_lethal)

    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_ben_lethal)))

#3  log-normal
if 3 in models:
    def prior_condition_log(par):
        if par[2]!=0:
            return True
        return False

    model = 3
    prior = [[0,1,'uniform'],[-1,1,'uniform'],[0,1,'uniform'],ga_prior]

    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_log)))

#4  log-normal truncated
if 4 in models:
    def prior_condition_log(par):
        if par[2]!=0:
            return True
        return False

    model = 4
    prior = [[0,1,'uniform'],[-1,1,'uniform'],[0,1,'uniform'],ga_prior]

    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_log)))

#5  5-spikes
if 5 in models:
    def prior_condition_spikes(par):
        if (par[0]+par[1]+par[2]+par[3])<=1:
            if par[4]<par[5]:
                return True
        return False

    model = 5
    prior = [[0,1,'uniform'],[0,1,'uniform'],[0,1,'uniform'],[0,1,'uniform'],
             [0,1,'uniform'],[0,1,'uniform'],[1,2,'uniform'],ga_prior]
    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_spikes)))

#6 neutral
if 6 in models:
    model = 6
    prior = [ga_prior]
    par = np.vstack((par,simulation_par(model,prior,N,maxPar,prior_condition_none)))


filename = sys.argv[1]
np.save(filename,par)
