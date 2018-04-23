# -*- coding: utf-8 -*-
"""
simulation_functions.py
author:     Eva Bons
created on:     Wed Mar 30 12:03:37 2016
last edited on: Wed Mar 30 12:03:37 2016

contains all functions required to perform a simulation of sequence evolution
"""
import numpy as np
import random
import sys
import os

def diverseSeq(seq,diversityLevel=0.075):
    seq1 = np.copy(seq)
    N_change = int(diversityLevel*len(seq))
    toChange = random.sample(xrange(0,len(seq)),N_change )
    seq1[toChange] = np.random.choice(range(4),N_change)
    return seq1

def createSeq(N):
    """creates a random sequence of N bases long with a base distrubtion of HIV """
    baseDistr = np.array([0.373,0.241,0.193,0.193]) #base distribution AGTC
    probDistr = np.cumsum(baseDistr) #take cumulative sum for indexing with random nrs
    seq = []
    for i in range(N):
        #find the base to use based on a random nr
        baseNr = random.random()
        baseIdx = np.where(probDistr>baseNr)[0][0]
        seq.append(baseIdx)
    return np.array(seq)

def nrMutations(gen,seq):
    """finds the number of unique mutations in (sample of) generation gen"""
    mut = []
    for i in gen:
        if i[0][0] != -1:
            for j in range(len(i[0])):
                mut.append(str(seq[i[0][j]])+'-'+str(i[0][j])+'-'+str(i[1][j]))
    return len(np.unique(mut))

def createFitnessTable(seq, model, par):
    '''creates a table with random fitness values according to the model provided:

        * 'neutral': no fintess effects

        * 'lethals only': only lethal and neutral mutations possible, par = [fl]

        * 'beneficials only': only beneficial (exponentially distributed and
           neutral mutations possible, par = [fb, beta]

        * 'lethals and beneficials': par = [fl,fb,beta]

        * 'spikes': several point-masses along MFED, par = [fl, f1, f2 ,..., b1, b2,...]

        * 'lognormal': lethals, rest lognormally distributed, par = [fl, mu,sigma]


        * 'lognormal truncated': lethals, rest lognormally distributed, no beneficials, par = [fl, mu,sigma]
    '''
    N = len(seq)

    fitness = np.zeros((4,N)) #create an array filled with 0 (all lethal)
    for (i, base) in zip(range(N),seq): #make sure the fitness benefit of the initial sequence is 1
            fitness[base,i] = 1
    toFill = np.where(fitness == 0)

    if model == 'neutral': #neutral model
        fitness = np.ones((4,N))
    elif model in ['lethals only','beneficials only','lethals and beneficials']:
        if model == 'lethals only':
            fl = par[0]
            fb = 0
            beta = 0
        elif model == 'beneficials only':
            fl = 0
            fb = par[0]
            beta = par[1]
        elif model == 'lethals and beneficials':
            fl = par[0]
            fb = par[1]
            beta = par[2]

        for i,j in zip(toFill[0],toFill[1]):
            randnr = random.random()
            if randnr > fl:
                if randnr > fl+fb: #neutral mutation
                    fitness[i,j] = 1
                else: # beneficial
                    fitness[i,j]=1+np.random.exponential(beta)
    elif model == 'spikes': #spikes
        fl = par[0]
        n_spikes = (len(par)-1)/2
        f_spikes = par[1:n_spikes+1]

        b_spikes = par[-n_spikes:]
        for i,j in zip(toFill[0],toFill[1]):
            randnr = random.random()
            if randnr > fl:
                if randnr > fl+sum(f_spikes): #neutral
                    fitness[i,j] = 1
                else:
                    prob = fl
                    for spike in range(n_spikes):
                        prob += f_spikes[spike]
                        if randnr < prob:
                            fitness[i,j] = b_spikes[spike]
                            break
    elif model in ['lognormal','lognormal truncated']:
        fl = par[0]
        mu = par[1]
        sigma = par[2]
        for i,j in zip(toFill[0],toFill[1]):
            randnr = random.random()
            if randnr > fl:
                fitness[i,j] = np.random.lognormal(mu,sigma)
                if model == 'lognormal truncated':
                    if fitness[i,j] > 1:
                        fitness[i,j] = 1
    elif model == '2exp':
        fd = par[0]
        fb = par[1]
        fl = 1-fd-fb
        ld = par[2]
        lb = par[3]
        for i,j in zip(toFill[0],toFill[1]):
            randnr = random.random()
            if randnr > fl:
                if randnr < fl + fd: #deleterious
                    fitness[i,j] = 1-np.random.exponential(1-ld)
                    if fitness[i,j] < 0:
                        fitness[i,j] = 0
                else:
                    fitness[i,j] = 1+np.random.exponential(lb)

    else:
        print 'unknown model'
        raise

    return fitness

def simulate(seq,fitnessTable,nrMut,nrSamples,initR0=6,GAincrease = 1,
             mutRate=2.16e-5,printGen=False):
    capacity = 10000 #carrying capacity of the populations: 10000
    mutRate =  (0.82*mutRate)+(0.18*GAincrease*mutRate) #adapt mutation rate for APOBEC increase in G-to-A mutations

    currentGen = [[[-1],[0]]] #current generation is the unchanged sequence
    #generations are stored as list of mutations, first list contains the location,
    #second list the base the original position mutated to (A=0,G=1,C=2,T=3)
    #e.g. two individuals with a mutation A->G at position 20 will be stored as [[[20],[1]],[[20],[1]]]

    nrGen = 50 #the maximum number of generations, to prevent an infinite loop
    nrMutOld = 0 #the number of mutations found in the previous generation, to find best match to data
    oldSample = currentGen
    for gen in range(nrGen):
        offsprings = [] #will contain the nr offspring of each current sequence
        newGenTemp = [] #will contain the indices of parents for next generation
        newGen = []
        sizeGeneration = len(currentGen)
        for i in range(sizeGeneration): #for all sequences in the current generation
            offsprings.append(getNrOffspring(currentGen[i],fitnessTable,initR0))
        if sum(offsprings) <= capacity:
            for i in range(len(offsprings)):
                newGenTemp = newGenTemp + [i for j in range(offsprings[i])]
        else:
            newGenTemp = np.random.choice(range(sizeGeneration),capacity,True,offsprings/np.sum(offsprings,dtype=float))


        if len(newGenTemp)>capacity:
            newGenTemp = random.sample(newGenTemp, capacity)

        for j in newGenTemp: #create (mutated) sequences for the 'surviving' offspring
            newGen.append(mutateSeq(seq,mutRate,currentGen[j],GAincrease))

        if len(newGen)<1: #if the population died out (highly unlikely after the first generation), redo the simulation
            return simulate(seq,fitnessTable,nrMut,nrSamples,initR0,GAincrease,
                            mutRate)

        if len(newGen)<nrSamples: #draw a sample for analysis
            currentSample = newGen
        else:
            currentSample = random.sample(newGen,nrSamples)
        nrMutCurrent = nrMutations(currentSample,seq)
        if nrMutCurrent >= nrMut: #once the desired nr mutations has been reached, end simulation
            if abs(nrMutCurrent-nrMut) < abs(nrMutOld-nrMut):
                return [currentSample, gen]
            else:
                return [oldSample,gen-1]
        currentGen = newGen
        oldSample = currentSample
        nrMutOld = nrMutCurrent
    if printGen:
        print gen
    if len(currentGen)>=nrSamples:
        return [random.sample(currentGen,nrSamples),gen]
    else:
        return [currentGen,gen]

def getNrOffspring(changes,fitnessTable,initR0=6):
    """returns the number of offspring of a sequence according to the fitness
    of that sequence"""
    fitness = 1
    if changes[0][0] != -1:
        for pos,base in zip(changes[0],changes[1]):
            fitness *= (fitnessTable[base,pos])
    try:
        return np.random.poisson(fitness*initR0)
    except ValueError:
        return fitness*initR0


def mutateSeq(seq,mutRate,mutations=[[-1],[0]],GAincrease = 1):
    """mutates a sequence (with existing mutations) of length N, according to
    the per base mutation rate"""
    baseSubst =  np.array([[ 0.54505451,  0.79257926,  0.90629063,  1. ],
                           [ 0.38363836,  0.81938194,  0.90849085,  1. ],
                           [ 0.21922192,  0.330033,    0.75747575,  1. ],
                           [ 0.18228177,  0.29717028,  0.54194581,  1. ]])
    N = len(seq)
    nrMut = np.random.binomial(N,mutRate) #draw how many mutations will take place
    successMut = 0
    while successMut < nrMut: #do the mutations one by one
        where = random.randrange(0,N) #draw where the mutation will take place
        randNr = random.random() #draw a random nr for base substitution
        if mutations[0][0]!=-1 and where in mutations[0]: #if the selected base has been mutated before
            where_mut = np.where(mutations[0]==where)[0][0] #find the position of the previous record of change
            toMutate = mutations[1][where_mut] #get the base to mutate
            toCheck= baseSubst[toMutate,:] #get the cumulative distribution of base substitutions
            mutatedIdx = np.where(toCheck>randNr)[0][0] #find the new base
            if (toMutate != mutatedIdx):
                if (toMutate == 1 and mutatedIdx == 0) and (GAincrease > 1 or random.random()<GAincrease): #accept GA:
                    mutations[1][where_mut] = mutatedIdx #replace with the new base to the mutations list
                    successMut+=1
                elif (toMutate != 1 or mutatedIdx != 0) and (GAincrease < 1 or random.random()<(1.0/GAincrease)): #accept other mutations
                    mutations[1][where_mut] = mutatedIdx #replace with the new base to the mutations list
                    successMut+=1
        else:
            toMutate = seq[where] #get the base to mutate
            toCheck= baseSubst[toMutate,:]#get the cumulative distribution of base substitutions
            mutatedIdx = np.where(toCheck>randNr)[0][0]#find the new base
            if (toMutate != mutatedIdx):
                if (toMutate == 1 and mutatedIdx == 0) and (GAincrease > 1 or random.random()<GAincrease): #accept GA
                    successMut+=1
                    if mutations[0][0]==-1: #if there has been no mutation in the sequence yet
                        mutations = np.array([[where],[mutatedIdx]]) #create an array with the new base substitution
                    else:
                        mutations = np.hstack((mutations, [[where],[mutatedIdx]])) #add base substitution to mutations list
                elif (toMutate != 1 or mutatedIdx != 0) and (GAincrease < 1 or random.random()<(1.0/GAincrease)): #accept other mutations
                    successMut+=1
                    if mutations[0][0]==-1: #if there has been no mutation in the sequence yet
                        mutations = np.array([[where],[mutatedIdx]]) #create an array with the new base substitution
                    else:
                        mutations = np.hstack((mutations, [[where],[mutatedIdx]])) #add base substitution to mutations list

    return mutations

def printIndividualStats(ID, results,f,nrSamples):
    unmutated = 0
    sites = []
    nrSamples = len(results)
    for seq in results:
        if seq[0][0]==-1:
            unmutated += 1
        else:
            for pos in seq[0]:
                if pos not in sites:
                    sites.append(pos)
    if unmutated != nrSamples:
        mean_selection = ( (len(sites)/2600.0) / ((nrSamples-unmutated)/float(nrSamples)) )
    else:
        mean_selection = 0
    mean_unmutated = unmutated/float(nrSamples)
    f.write('#{}\t{}\t{}\t{}\n'.format(ID,mean_unmutated,mean_selection,nrSamples))

def printMutations(gen,seq,sim,f):
    """prints all mutation occuring in the generations in tab-seperated format"""
    for (i,nr) in zip(gen,range(len(gen))):
        if i[0][0] != -1:
            for j in range(len(i[0])):
                if i[1][j] != seq[i[0][j]]:
                    f.write(str(seq[i[0][j]])+'-'+str(i[0][j]) + '-' +
                               str(i[1][j])+'\t'+str(sim)+'\t'+str(nr)+'\n')

def simulate_and_print(seq,fitnessTable,nrMut,nrSamples,ID,fileDirect,printGen=False,GAincrease=1):
    f = open(fileDirect+'/out_'+str(ID),'w')
    results = simulate(seq,fitnessTable,nrMut,nrSamples,printGen=printGen,GAincrease=GAincrease)
    printIndividualStats(ID,results[0],f,nrSamples)
    printMutations(results[0],seq,ID,f)
    f.write('\n')
    f.close()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    seq = np.zeros(1000)
    print seq
    fitness = createFitnessTable(seq, model = '2exp', par = [0.7,0.1,0.8,0.05])
    plt.hist(fitness[1,:])
    plt.show()
