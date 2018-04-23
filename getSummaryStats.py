"""
getSummaryStats.py
author:     Eva Bons
created on:     Feb 25 2016
last edited on: Feb 25 2016

usage: python getSummaryStats.py infile [>outfile]

from a list of mutations create by simulate.py, get the summary statistics
necessary for ABC_MCMC
prints:
    mean of fraction shared per individual
    median of fraction shared per individual
    mean-median of fraction shared per individual
    variance of fraction shared per individual
    mutations occuring in 1 individual
    mutations occuring in 2 individuals
    mutations occuring in 3-5 individuals
    mutations occuring in 6-10 individuals
    mutations occuring in >10 individuals
    area under cumulative fraction shared curve
    mean fraction of unmutated seqeunces in the samples
    mean selection coefficient (fraction mutated sites/fraction mutated seq)
    mean Poisson estimate of intersequence HD
    fraction of mutations that is G-to-A
"""


import numpy as np
import sys

def listMut(results,seq):
    all_mutations = []
    for inidividual in range(len(results)):
        all_mutations.append([])
        mutations = results[str(inidividual)][0]
        for s in mutations:
            if s[0][0]!=-1:
                for mut,i in zip(s[0],range(len(s[0]))):
                    if seq[mut] != s[1][i]:
                        all_mutations[-1].append(str(seq[mut])+'-'+str(mut)+'-'+str(s[1][i]))
    return all_mutations

def shared_per_ind(results,seq=0,asList=False):
    if not asList:
        all_mutations = listMut(results,seq)
    else:
        all_mutations = results
    fractIndnrShared = []
    for i in range(len(all_mutations)):
        uniqueMut = np.unique(all_mutations[i])
        nrShared = 0
        for j in uniqueMut: #go through all unique mutations in this simulation
            shared = False
            for other in (all_mutations[:i]+ all_mutations[i+1:]): #go through all other individuals
                if j in other:
                    if not shared: #only add each shared mutation once
                        nrShared += 1.0
                        shared = True
        if len(uniqueMut)>0:
            fractIndnrShared.append(nrShared/len(uniqueMut))
        else:
            fractIndnrShared.append(0)
    return fractIndnrShared

def nr_occuring(results, seq='',asList=False):
    if not asList:
        listShared = listMut(results,seq)
    else:
        listShared = results
    muts = []
    nrs = []
    for i in listShared:
        for j in np.unique(i):
            if j in muts:
                nrs[muts.index(j)]+=1
            else:
                muts.append(j)
                nrs.append(1)
    counts = np.histogram(nrs,bins=[0,1.1,2.1,5.1,10.1,max(100,max(nrs)+1)])
    return counts[0]

def area_under_curve(results, seq='',asList=False):
    listShared = np.cumsum(sorted(shared_per_ind(results,seq,asList)))
    return np.trapz(listShared,dx=1.0/len(listShared))

def area_under_curve2(results,seq='',asList=True):
    listShared=np.sort(shared_per_ind(results,seq,asList))
    previous=0
    counts=0
    total = 0
    for i in listShared:
        if i == previous:
            counts+=1
        else:
            total += ((float(counts)/len(listShared))*(i-previous))
            previous = i
    return total




def simple_stats(results, seq='',asList=False):
    listShared = shared_per_ind(results,seq,asList)
    stats = {}
    stats['mean'] = np.mean(listShared)
    stats['median'] = np.median(listShared)
    stats['variance'] = np.var(listShared)
    return stats


def individual_stats(results):
    mean_unmutated = 0
    mean_selection = 0
    for individual in range(len(results)):
        ind_mut = results[str(individual)][0]
        unmutated = 0
        sites = []
        for seq in ind_mut:
            if seq[0][0]==-1:
                unmutated += 1
            else:
                for pos in seq[0]:
                    if pos not in sites:
                        sites.append(pos)
        mean_selection += ( (len(sites)/2600.0) / ((len(ind_mut)-unmutated)/float(len(ind_mut))) ) / len(results)
        mean_unmutated += (unmutated/float(len(ind_mut)))/(len(results))
    return(mean_unmutated,mean_selection)

def hammingDist(seq1,seq2):
    HD = 0

    posA=[]
    baseA=[]
    for i in seq1:
        A = i.split('-')
        posA.append(A[1])
        baseA.append(A[2])

    posB=[]
    baseB=[]
    for i in seq2:
        B = i.split('-')
        posB.append(B[1])
        baseB.append(B[2])
    for [pos,base] in zip(posA,baseA):
        if pos in posB:
           if np.array(baseB)[np.where(np.array(posB)==pos)]!=base:
               HD += 0.5
        else:
            HD+=1
    return HD

def meanHD(gen,seq,nrSeq):
    distances = HDDistribution(gen,seq,nrSeq)
    return np.mean(distances)

def HDDistribution(gen,seq,nrSamples):
    size = len(gen)
    distances = []
    for i in range(nrSamples):
        seqA = np.array(gen)[np.where(np.array(seq)==i)]
        for j in range(i+1,nrSamples):
            seqB = np.array(gen)[np.where(np.array(seq)==j)]
            HD = hammingDist(seqA,seqB)
            distances.append(HD)
    return distances

def GtoAfraction(muts):
    GtoA = 0.0
    total = 0
    for patient in muts:
        for mutation in patient:
            if mutation[0] in ['G','1'] and mutation[-1] in ['A','0']:
                GtoA+=1
            total+=1
    return GtoA/total

def get_stats(arg):
    f = open(arg)

    ind = ''
    nrSamples = []
    mut= []
    HD = []
    frac_unmut = []
    selection = []
    for line in f.readlines():
        if line != '\n':
            fields = line.split()
            if '#' in line:
                frac_unmut.append(float(fields[1]))
                selection.append(float(fields[2]))
                nrSamples.append(int(fields[3]))
            else:
                if fields[1]!=ind:
                    if len(mut)>0:
                        HD.append(meanHD(mut[-1],seq,nrSamples[-1]))
                    mut.append([fields[0]])
                    ind=fields[1]
                    seq = [int(fields[2])]
                else:
                    mut[-1].append(fields[0])
                    seq.append(int(fields[2]))
    f.close()
    stats = []
    simple_stat = simple_stats(mut,asList=True)
    stats.append(simple_stat['mean'])
    stats.append(simple_stat['median'])
    stats.append(simple_stat['mean']-simple_stat['median'])
    stats.append(simple_stat['variance'])

    counts = nr_occuring(mut,asList=True)
    stats.append(counts[0]) #singletons
    stats.append(counts[1]) #occuring twice
    stats.append(counts[2]) #occurring 2-5 times
    stats.append(counts[3]) #occuring 5-10 times
    stats.append(counts[4]) #occuring >10 times

    stats.append(area_under_curve2(mut,asList=True))

    stats.append(np.mean(frac_unmut))
    stats.append(np.mean(selection))
    stats.append(np.nanmean(HD))
    stats.append(GtoAfraction(mut))
    return stats

if __name__ == '__main__':
    stats = main(sys.argv[1])
    for i in stats:
        print i
