#!/usr/bin/env python2

import hifive
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

#generating the heatmap
hic=hifive.HiC('hiproject', 'r')

data = hic.cis_heatmap('chr13', 1000000, datatype='fend', arraytype='full', diagonalincluded=True)

zerofilter = np.where(data[:,:,0:2] > 0)

nozero_data = data[zerofilter[0], zerofilter[1], 0]
nozero_data2 = data[zerofilter[0], zerofilter[1], 1]

enrichment = nozero_data / nozero_data2

heatmap_array = np.zeros((1193,1193))

for i in range(0, len(enrichment)):
    heatmap_array[zerofilter[0][i]][zerofilter[1][i]] = np.log(enrichment[i])



fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap="Blues", aspect = 'auto', interpolation='none')#the basic heatmap

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel('log Expression', fontsize = 12, rotation=-90, va="bottom")
ax.set_title("Chromosome 13 Heatmap", fontsize = 15)

fig.set_size_inches(10.5, 10.5)
fig.tight_layout()

plt.savefig('Heatmap_week9.png')

#compartment analysis
Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
Comp.write_eigen_scores('hic_comp.bed')

fig, ax = plt.subplots()
ax.plot(Comp.positions['chr13'],Comp.eigenv['chr13'])
ax.set_title("Compartment plot", fontsize = 15)
plt.xlabel('Position on chromosome13', fontsize = 14)
plt.ylabel('Eigenv', fontsize = 14)
plt.savefig('Compartment_plot.png')

#Violin plots for expression
pos = np.genfromtxt("pos_expression.bed", usecols=(4))
neg = np.genfromtxt("neg_expression.bed", usecols=(4))

compartment_labels = ['B', 'A']
fig, ax = plt.subplots()
ax.violinplot(pos, positions=[0], showextrema=False)
ax.violinplot(neg, positions=[1], showextrema=False)

ax.set_title("Violin Plot for Compartments", fontsize = 15)
ax.set_xticks([0,1])
ax.set_xticklabels(['B', 'A'])
ax.set_yscale('log')
plt.xlabel('Compartments', fontsize = 14)
plt.ylabel('Gene Expression', fontsize = 14)
plt.savefig('Violin_plot.png')

#Repression analysis
positive_rep = []
negative_rep = []

bw = pyBigWig.open('data/WT_H3K27me3.bw')
fs_neg = open('neg_expression.bed', 'r')
fs_pos = open('pos_expression.bed', 'r')

for line in fs_neg:
    #split each by whitespace into columns
    line_split = line.split()
    negative_rep.append(bw.stats('chr13', int(line_split[1]), int(line_split[2]), type='sum'))
    
for line in fs_pos:
    #split each by whitespace into columns
    line_split = line.split()
    positive_rep.append(bw.stats('chr13', int(line_split[1]), int(line_split[2]), type='sum'))

pos_flat = []
neg_flat = []
exp_pos = pos.tolist()
exp_neg = neg.tolist()

for sublist in positive_rep:
    for i in sublist:
        pos_flat.append(i)

for sublist in negative_rep:
    for i in sublist:
        neg_flat.append(i)
        
[0 if value is None else value for value in pos_flat]
[0 if value is None else value for value in neg_flat]

fig, ax = plt.subplots(ncols=2, figsize = (10,5))
ax[0].scatter(exp_pos,pos_flat, alpha = 0.5)
ax[1].scatter(exp_neg,neg_flat, alpha = 0.5)
ax[0].set_title("Expression vs Repression \n in B compartment", fontsize = 15)
ax[1].set_title("Expression vs Repression \n in A compartment", fontsize = 15)
ax[0].set_ylabel('Repression')
ax[1].set_ylabel('Repression')
ax[0].set_xlabel('Expression')
ax[1].set_xlabel('Expression')

plt.tight_layout()
plt.savefig('REPGraph.png')
