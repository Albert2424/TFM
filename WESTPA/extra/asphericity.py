#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 12:56:12 2023

@author: albert

asphericity
"""

import os
import shutil
import mdtraj as md
import re
import numpy as np
from cluster import *
from analyse import *
from MDAnalysis import transformations
import MDAnalysis
from MDAnalysis.tests.datafiles import PDB

def get_clusts(clusters,fasta,prot,frame):
    
    # Get the centers of each cluster
    biggest = 0
    for clust in clusters[frame]:
        size = clusters[frame][clust]['size']
        # print(size)
        if size > biggest:
            chain_per_clust = clusters[frame][clust]['chains'] #chains of the clust 
            positions = np.array(pos[frame])[chain_per_clust] #CM of the clust

            biggest = size
        
    return chain_per_clust

def config_list(windows,n_chains,seq):
    """
    Generates a list of all the configurations in the /config directory.
    -------------------------------------------------------------
    
    INPUT:
    -------
        windows: integer. Number of initial configurations that the /config 
                 directory contains.             
        
    OUTPUT:
    --------
        config: list. List of configurations that will be graphed.
        
        expected: list. list containing the expected values of the cluster sizes.

    """
    # config = []
    # expected = []
    # for c in range(windows):
    #     perc = f'{int(int(c/(windows-1)*100)):02d}' #name of the end of the config
    #     file = perc #initial configurations
    #     config.append(file)
    config = []
    with open('bstates_'+seq+'/bstates.txt','r') as f:
        for line in f:
            config.append(line.split()[-1])
    return config


if __name__ == '__main__':
    
    seq = 'WT'   
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    n_chains = 150
    L = 343
    N = len(fasta)
    
    config = config_list(10,n_chains,seq)
    # config=['00','11','22','33','44','56','67','78','89','100']
    
    for i in config:
        filename = f'bstates_{seq}/{i}/top.pdb'
    
        t = md.load(filename)
        
        ipos = get_initial_pos(filename)
        prot = protein(ipos,n_chains)
        pos = get_points(fasta,prot)
        cl = clust(pos,17.,L,2,fasta,prot)#get the clusters
        
        chains = get_clusts(cl,fasta,prot,'frame 0') #get clust chains
        chains = np.array(chains)
        cluster_resi = t.top.select("chainid "+' '.join([str(j) for j in chains]))
        t = t.atom_slice(cluster_resi)
        bonds=np.array([(i,i+1) for i in np.arange(t.n_atoms-1)], dtype=np.int32)
        t.make_molecules_whole(inplace=True,sorted_bonds=bonds)
        t.center_coordinates()
      
        

        print(f'Cluster with {len(chains)} chains and realtive shape anisotropy: {md.relative_shape_antisotropy(t)[0]:.3f}')
        
#%%

import matplotlib.pyplot as plt
import matplotlib.cm as cm

c_wt250 = [5,26,55,83,110,137,167,194,221,250]
c_wt150 = [3,17,32,48,65,84,100,114,132,148]
k250_wt = [0.440,0.021,0.074,0.036,0.095,0.073,0.078,0.152,0.094,0.093]
k150_wt = [0.610,0.044,0.067,0.038,0.035,0.102,0.071,0.058,0.118,0.105]


c_s250 = [4,26,54,81,109,137,166,194,222,250]
c_s150 = [3,15,32,49,67,81,99,114,132,148]
k250_s = [0.590,0.075,0.134,0.196,0.098,0.134,0.109,0.112,0.128,0.136]
k150_s = [0.248,0.191,0.072,0.050,0.113,0.018,0.230,0.146,0.1,0.087]

col = cm.viridis(np.linspace(0, 1, 4))

plt.subplots(2,1,figsize=(8,7))

plt.subplot(211)
plt.plot(c_wt250[2:],k250_wt[2:], label = r'WT ($c_2$)',color=col[0],marker='o',
markersize=5,markeredgewidth=0.5,markeredgecolor='black')
plt.plot(c_s250[2:],k250_s[2:], label = r'Shuffle ($c_2$)',color=col[1],marker='o',
markersize=5,markeredgewidth=0.5,markeredgecolor='black')
plt.xlabel('CS',fontsize=20)
plt.ylabel(r'$\kappa^2$ (nm)',fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.ylim(0.0,0.3)
plt.legend(fontsize=20)

plt.subplot(212)
plt.plot(c_wt150[2:],k150_wt[2:], label = r'WT ($c_1$)',color=col[2],marker='o',
markersize=5,markeredgewidth=0.5,markeredgecolor='black')
plt.plot(c_s150[2:],k150_s[2:], label = r'Shuffle ($c_1$)',color=col[3],marker='o',
markersize=5,markeredgewidth=0.5,markeredgecolor='black')


plt.xlabel('CS',fontsize=20)
plt.ylabel(r'$\kappa^2$ (nm)',fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.legend(fontsize=20)
plt.tight_layout()

plt.savefig('rsa.pdf')
plt.show()