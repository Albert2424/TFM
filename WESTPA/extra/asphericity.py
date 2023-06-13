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


if __name__ == '__main__':
    
    filename = 'conf/100/top.pdb'
    seq = 'WT'
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    n_chains = 75
    L = 273
    N = len(fasta)
    
    t = md.load(filename)
    
    ipos = get_initial_pos(filename)
    prot = protein(ipos,n_chains)
    pos = get_points(fasta,prot)
    cl = clust(pos,10.,L,2,fasta,prot)#get the clusters
    
    chains = get_clusts(cl,fasta,prot,'frame 0') #get clust chains
    chains = np.array(chains)
    # clust_chains = []
    # for i in chains: #get clust atoms 
    #     i += 1
    #     # print(ipos[0][(i-1)*N:(i)*N])
    #     clust_chains += ipos[0][(i-1)*N:(i)*N] #select chain
    
    t = t.atom_slice(chains)
    print(f'Cluster with {len(chains)} chains and asphericity: {md.asphericity(t)[0]}')
        
