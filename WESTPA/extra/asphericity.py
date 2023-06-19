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

def config_list(windows,n_chains):
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
    config = []
    expected = []
    for c in range(windows):
        perc = f'{int(int(c/(windows-1)*100)):02d}' #name of the end of the config
        file = perc #initial configurations
        config.append(file)
    return config


if __name__ == '__main__':
    
    seq = 'WT'   
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    n_chains = 150
    L = 343
    N = len(fasta)
    
    config = config_list(10,n_chains)
    # config=['00','11','22','33','44','56','67','78','89','100']
    
    for i in config:
        filename = f'bstates_WT/{i}/top.pdb'
    
        t = md.load(filename)
        
        ipos = get_initial_pos(filename)
        prot = protein(ipos,n_chains)
        pos = get_points(fasta,prot)
        cl = clust(pos,17.,L,2,fasta,prot)#get the clusters
        
        chains = get_clusts(cl,fasta,prot,'frame 0') #get clust chains
        chains = np.array(chains)
        cluster_resi = t.top.select("chainid "+' '.join([str(i) for i in chains]))
        t = t.atom_slice(cluster_resi)
        bonds=np.array([(i,i+1) for i in np.arange(t.n_atoms-1)], dtype=np.int32)
        t.make_molecules_whole(inplace=True,sorted_bonds=bonds)
        t.center_coordinates()
      
        

        print(f'Cluster with {len(chains)} chains and realtive shape anisotropy: {md.relative_shape_antisotropy(t)[0]:.3f}')
        
