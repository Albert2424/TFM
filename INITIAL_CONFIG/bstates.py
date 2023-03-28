#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 12:01:53 2023

@author: albert
"""

"""CLUSTER DETECTION"""

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN,KMeans
from analyse import *
from argparse import ArgumentParser
import os
import shutil

# parser = ArgumentParser()
# parser.add_argument('--name',nargs='?',const='', type=str)
# args = parser.parse_args()

name='WT'


def get_traj(filename,top):
    
    """
    Provides the position of all the particles of the system in the
    configuration provided by the filename on every frame of the dcd.
    -------------------------------------------------------------
    
    INPUT:
    -------
        filename: .dcd file containing the frames of the simulation.
        
        top: .pdb file containing the configuration.
        
        
    OUTPUT:
    --------
        traj: array of float of size (frames, npar, 3). Returns a 3D array where
              the first dimesnion is the number of frames, the second corresponds
              to the number of particles of the system and the third is the 
              xyz coordinates of each particle.                  
    
    """
    
    file=md.load_dcd(filename,top)
    traj=file.xyz
    return traj

def get_initial_pos(filename):
    
    """
    Provides the position of all the particles of the system in the
    configuration provided by the filename.
    -------------------------------------------------------------
    
    INPUT:
    -------
        filename: .pdb file containing the configuration.
        
    OUTPUT:
    --------
        ipos: array of float of size (1, npar, 3). Returns a 3D array where
              the first dimesnion is 1 (only 1 frame), the second corresponds
              to the number of particles of the system and the third is the 
              xyz coordinates of each particle.                  
    
    """
    
    file=md.load_pdb(filename)
    ipos=file.xyz
    return ipos

def protein(traj,n_chains):
    
    """
    
    Divides the system in chains for each frame.
    -------------------------------------------------------------
    
    INPUT:
    -------
        traj: array of float of size (frames, npar, 3). Returns a 3D array where
              the first dimesnion is the number of frames, the second corresponds
              to the number of particles of the system and the third is the 
              xyz coordinates of each particle.

        n_chains: integer. Number of chains of the system.                  
    
        
    OUTPUT:
    --------
        prot: dictionary containing all frames of the trajectory.                 
    
    """
    
    nprot=len(traj[0,:,:])
    fasta=nprot/n_chains
    
    prot={}
    
    for frame in range(len(traj)):
        prot['frame {:}'.format(frame)] = []
        for i in range(n_chains):
            prot['frame {:}'.format(frame)].append(traj[frame,int(fasta*(i)):int(fasta*(i+1)),:])
        
    return prot

def CM(fasta,prot):
    
    """
    
    Computes the center of mass of each chain for every frame. The
    mass of each aminoacid is in Daltons (1/N_A).
    -------------------------------------------------------------
    
    INPUT:
    -------
        fasta: integer. Number of particles per chain.               
        
        prot: dictionary containing all frames of the trajectory.
   
    OUTPUT:
    --------
        r_CM: dictionary containing  for all frames the center of mass of 
        each chain.                 
    
    """
    
    aa_mass = {
    'A': 71.03711, 'G': 57.02146, 'M': 131.04049, 'S': 87.03203,
    'C': 103.00919, 'H': 137.05891, 'N': 114.04293, 'T': 101.04768,
    'D': 115.02694, 'I': 113.08406, 'P': 97.05276, 'V': 99.06841,
    'E': 129.04259, 'K': 128.09496, 'Q': 128.05858, 'W': 186.07931,
    'F': 147.06841, 'L': 113.08406, 'R': 156.10111, 'Y': 163.06333
    } #mass in Dalton (1/N_A)
    
    mass=[]
    
    for amin in fasta:
        mass.append(aa_mass[amin])
        
    M = np.sum(mass)
    
    r_CM = {}
    for frame in prot:
        r_CM[frame] = []
        for p in prot[frame]:
            r_CM[frame].append(np.dot(np.array(mass), p) / M)
        
    return r_CM

def wolf_algorithm(X, radius, min_samples=2):
    
    """
    Uses the wolf algorithm to sample the different clusters of the system.
    -------------------------------------------------------------
    
    INPUT:
    -------
        X: array (npar,3). Array containing the position of every center of mass
                           of the frame.                
        
        radius: float. Indicate the radius at which a particle is part of a cluster.
        
        min_samples: integer. Minimum number of particles that hte algorithm considers
                     as a cluster (default is 2).
   
    OUTPUT:
    --------
        dbscan.labels_: 1D array containing the label of each particle. Particles
                        with the same label correspond to the same cluster. Particles
                        with -1 label are not part of any cluster.
    
    """
    
    dbscan = DBSCAN(eps=radius, min_samples=min_samples)
    dbscan.fit(X)
    return dbscan.labels_


def clust(pos,dist,min_size):
    
    """
    Detects clusters formed in every frame of the simulation.
    -------------------------------------------------------------
    
    INPUT:
    -------
        pos: dict with keys of the kind 'frame 0', 'frame 1',... dictionary 
        containing all frames of the trajectory.              
        
        dist: float. Indicate the radius at which a particle is part of a cluster.
        
        min_size: integer. Minimum number of particles that hte algorithm considers
                     as a cluster (default is 2).
   
    OUTPUT:
    --------
        clusters: dictionary containing every cluster of the frame. Each cluster 
        has its own dictionary that contains the chains that are part of it, the
        position of the center of mass of said chains and the size and radius
        of the cluster:
            
            clusters = {0:{'chains':[0,1,2...],
                           'pos':[np.array(3),np.array(3),..],
                           'size': 14,
                           'radius': 1.453},
                        1:{...}
                        }
            
        centers: array containing the xyz coordinates of every cluster center.
                        
    
    """
    
    clusters={}
    
    for frame in pos:
        clusters[frame] = {}
        
        #calculate the different clusters
        labels = wolf_algorithm(pos[frame], dist, min_size)
        
        #add the clusters into the dictionary
        for label in range(len(labels)):
            
            if labels[label] != -1 :
                if labels[label] in clusters[frame]:
                    clusters[frame][labels[label]]['chains'].append(label)
                    clusters[frame][labels[label]]['pos'].append(np.array(pos[frame])[label])
                    
                    clusters[frame][labels[label]]['size']+=1
                else: 
                    clusters[frame][labels[label]]={'chains':[label],'pos':[np.array(pos[frame])[label]],'size':1,'rad':0}
                    
        
        #RADIUS OF THE CLUSTER
        
        # Get the centers of each cluster
        centers = []
        for frame in clusters:
            for clust in clusters[frame]:
                centers.append([np.mean(np.array(clusters[frame][clust]['pos'])[:,0]),
                                np.mean(np.array(clusters[frame][clust]['pos'])[:,1]),
                                np.mean(np.array(clusters[frame][clust]['pos'])[:,2])])
                  
                distances = np.linalg.norm(np.array(clusters[frame][clust]['pos']) - centers[clust], axis=1)
                clusters[frame][clust]['rad'] = np.max(distances)#add the radius        
                    
        return clusters,np.array(centers)

def directories(windows,name):
    
    print('creating /bstates directory')

    #clean if there was a previous directory
    
    try:
        shutil.rmtree('bstates')
    except OSError:
        pass
    
    #create the bstates directory
    
    os.mkdir('bstates')
    
    #initialize the protein
    
    proteins = initProteins()
    fasta_WT = proteins.loc[name].fasta

    for config in range(windows):
        perc = f'{config/(windows-1):.2f}' #name of the end of the config
        file = 'top_eq_'+perc+'.pdb' #initial configurations
        directory = f'{int(float(perc)*100):02d}' #directories
        
        os.mkdir('bstates/'+directory)
        print('/'+directory)
        
        shutil.copyfile('config/'+file,'bstates/'+directory+'/'+file) #copy the file into the dir
        
        file = 'config/'+file

        ipos=get_initial_pos(file)

        prot=protein(ipos,100)

        pos=CM(fasta_WT,prot)

        cl,centers=clust(pos,7.,2)

        print('Number of clusters: ',len(cl['frame 0']))
        rad = []
        for i in cl['frame 0']:
            print(f'cluster {i:} size: {cl["frame 0"][i]["size"]:} and radius: {cl["frame 0"][i]["rad"]:.6f}')
            
            rad.append(cl["frame 0"][i]["rad"])   
        with open('bstates/'+directory+'/pcoord.ini','w') as f:
            try:
                f.write(str(np.max(rad)))
            except ValueError:
                f.write(str(0.))
            pass
                
        with open('bstates/pcoord.ini','a') as f:
            try:
                f.write(str(np.max(rad))+'\n')
            except ValueError:
                f.write(str(0.)+'\n')
            pass

        with open('bstates/bstates.txt','a') as f:
            f.write(str(config)+' '+str(round(1/windows,2))+' '+directory+'\n')
                
        print('')

        
        
        

        
        
        
#%%

windows = 6
name = 'WT'

directories(windows,name)

