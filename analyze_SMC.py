#imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scats
from collections import Counter

from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#settings
n_iter = 7
max_sims = [60,50,100,70,80,130,200]
max_dist = [5,2.2,1.3,0.8,0.7,0.6,0.5,0.4]
model_names = ['beneficials only', 'lethals only', 'lethals and beneficials',
               'lognormal', 'lognormal truncated', 'spikes', 'neutral']
par_names = [[r'$f_b$',r'$\beta$',r'G$\rightarrow$A'],[r'$f_l$',r'G$\rightarrow$A'],
             [r'$f_l$',r'$f_b$',r'$\beta$',r'G$\rightarrow$A'],[r'$f_l$',r'$\mu$',r'$\sigma$',r'G$\rightarrow$A'],
             [r'$f_l$',r'$\mu$',r'$\sigma$',r'G$\rightarrow$A'],
             [r'$f_l$',r'$f_{d_1}$',r'$f_{d_2}$',r'$f_b$',r'$b_{d_1}$',r'$b_{d_2}$',r'$b_b$',r'G$\rightarrow$A'],
             [r'G$\rightarrow$A']]

par_ranges = [[[0,1],[0,2],[-10,np.log(160)]],[[0,1],[-10,np.log(160)]],[[0,1],[0,1],[0,2],[-10,np.log(160)]],
              [[0,1],[-1,1],[0,1],[-10,np.log(160)]],[[0,1],[-1,1],[0,1],[-10,np.log(160)]],
              [[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[1,2],[-10,np.log(160)]],[[-10,np.log(160)]]]

colors= ['#FDEB2D','#0887B0','#AD1B28','#F638F2','#1604D6','#87D96A','#E98A7B']
markers = ['d','x','+','*','v','^','o' ]


#read in the data
stats = {}

for j in range(1,n_iter+1):
    first=True
    for i in range(1,max_sims[j-1]):
        if first:
            try:
                stats[j] = np.load('stats_set_{}/stats_{}.npy'.format(j,i))
                first=False
            except IOError:
                first=True
                print 'not found: {}'.format(i)
        else:
            try:
                stats[j] = np.vstack((stats[j], np.load('stats_set_{}/stats_{}.npy'.format(j,i))))
            except IOError:
                print 'not found: {}'.format(i)

    n_accepted = sum(stats[j][:,-1]<max_dist[j])

#SMC model probabilities
all_iterations = {i: [1.0/7] for i in range(7)}

plt.figure(figsize=[10,3])
for it in range(1,n_iter+1):
    accepted_models = stats[it][stats[it][:,-1]<max_dist[it],0]
    total_accepted = len(accepted_models)
    accepted_counts = Counter(accepted_models)

    for i in accepted_counts:
        all_iterations[i].append(1.0*accepted_counts[i]/total_accepted)

for i in all_iterations:
    all_iterations[i] = all_iterations[i] + [0 for j in range(1+n_iter-len(all_iterations[i]))]
    plt.scatter(range(len(all_iterations[i])), all_iterations[i],
                color=colors[int(i)],marker=markers[int(i)],s=100,
                label=model_names[i],alpha=0.5)
    plt.plot(range(len(all_iterations[i])), all_iterations[i], color=colors[int(i)],alpha=0.5)
    plt.text(7,all_iterations[i][-1],'{:.3f}'.format(all_iterations[i][-1]),
             ha='right',color=colors[int(i)])
ax=plt.axis()
plt.axis([0,n_iter,0,ax[-1]])

# Put a legend to the right of the current axis
lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Iteration')
plt.ylabel('Model probability')
plt.title('model probabilities per iteration')
plt.show()

#posterior distributions of last iteration, with best estimate and HPD intervals
for model in range(7):
    n_par = len(par_names[model])
    stats_here = stats[n_iter][(stats[n_iter][:,0]==model) & (stats[n_iter][:,-1]<max_dist[-1])]
    par_here = stats_here[:,1:n_par+1]
    par_here[:,-1] = np.log(par_here[:,-1])

    for i in range(n_par):
         #plot histogram
        single_par = par_here[:,i]
        if len(single_par) < 100:
            plot=False
        else:
            plot = True
        if plot:
            plt.subplot(1,n_par,i+1)

            bins = np.linspace(min(single_par), max(single_par),10)
            plt.hist(single_par,edgecolor='white',alpha=0.2,normed=True, color='grey',bins=bins)

            #calculate and plot density
            xs = np.linspace(min(single_par), max(single_par),1000)
            density = scats.gaussian_kde(single_par)
            ys = density.evaluate(xs)
            #calculate 95%HDP
            height = 0
            step = max(ys)/1000
            while  height < max(ys):
                if 1.0*sum(ys[ys>height])/sum(ys) < 0.95:
                    limits = (min(xs[ys>height]),max(xs[ys>height]))
                    break
                height+=step

            density_stats = {'xs': xs, 'ys': ys, 'HPD': limits, 'max': xs[np.argmax(ys)]}
            plt.plot(density_stats['xs'],density_stats['ys'],color='grey')

            #calculate and plot HPD intervals
            limits = density_stats['HPD']
            plt.axvspan(limits[0], limits[1], alpha=0.1, color='black')

            #calculate the maximum
            plt.axvline(density_stats['max'],ls= '--', color='grey')
            plt.text(density_stats['max'],max(density_stats['ys']),'{:.2f}'.format(density_stats['max']))

            plt.xlabel(par_names[model][i])
            plt.ylabel('density')
            plt.title('{:.3f} ({:.3f},{:.3f})'.format(density_stats['max'], limits[0],limits[1]))
            ax = plt.axis()
            plt.axis([par_ranges[model][i][0],par_ranges[model][i][1],ax[2],ax[3]])
    if plot:
        plt.suptitle(model_names[model])
        plt.show()
