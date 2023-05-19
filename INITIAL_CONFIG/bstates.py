#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 12:01:53 2023

@author: albert
"""

"""CLUSTER DETECTION"""

import numpy as np
import mdtraj as md
# import MDAnalysis
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN,KMeans
from analyse_ini import *
from argparse import ArgumentParser
import os
import shutil
from scipy.spatial import distance

parser = ArgumentParser()
parser.add_argument('--seq',nargs='?',const='', type=str)
parser.add_argument('--windows',nargs='?',const='', type=int)
parser.add_argument('--rc',nargs='?',const='', type=float)
parser.add_argument('--L',nargs='?',const='', type=float)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
args = parser.parse_args()


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
    ipos = file.xyz
    
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

def clust_detection(X, radius, L,min_samples=2):
    
    """
    Uses the DBSCAN algorithm to sample the different clusters of the system.
    -------------------------------------------------------------
    
    INPUT:
    -------
        X: array (npar,3). Array containing the position of every center of mass
                           of the frame.                
        
        radius: float. Indicate the radius at which a particle is part of a cluster.
        
        min_samples: integer. Minimum number of particles that hte algorithm considers
                     as a cluster (default is 2).
        L: float. Size of the simulation box.
   
    OUTPUT:
    --------
        dbscan.labels_: 1D array containing the label of each particle. Particles
                        with the same label correspond to the same cluster. Particles
                        with -1 label are not part of any cluster.
    
    """
    X=np.array(X) #avoids [array(),array(),...]
    def my_pdist(x,y,L):
        dx = np.abs(x-y)
        dx = np.abs(dx - np.rint(dx/L)*L)
        return np.linalg.norm(dx)


    dbscan = DBSCAN(eps=radius, metric=my_pdist, metric_params={'L':L}, min_samples=min_samples).fit(X)
    return dbscan.labels_


def clust(pos,dist,L,min_size):
    
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
        L: float. Size of the simulation box.
   
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
        labels = clust_detection(pos[frame], dist,L, min_size)
        
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
                  
                # distances = np.linalg.norm(np.array(clusters[frame][clust]['pos']) - centers[clust], axis=1)
                
                #rotation radius
                # clusters[frame][clust]['rad'] = np.max(distances)#add the radius 
                
                #hidrodynamic radius
                clusters[frame][clust]['rad'] = 1/(np.average(1/distance.pdist(np.array(clusters[frame][clust]['pos']))))
                    
        return clusters,np.array(centers)

def directories(windows,name,rc,L,n_chains):
    
    """
    Creates the bstates directory necessary for the performance
    of a westpa simulation. Notice that the pcoord of the system
    is the radius of the biguer cluster.
    -------------------------------------------------------------
    
    INPUT:
    -------
        windows: integer. Number of initial configurations that the /config 
                 directory contains.              
        
        name: str. Name of the protein used in the generation of the initial config.
        
        L: float. Size of the simulation box.
        
        n_chains: int. Number of chains of the simulation.
        
        rc: float. Minimum distance to consider a particle in the cluster.
        
       
    OUTPUT:
    --------
        /bstates: directory contining all the necessary files for the westpa 
                  simulation bstate directory (use the generated directory on 
                  your westpa instead of the default). 
                
    """
    
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
    
    bins=[]

    for config in range(windows):
        perc = f'{config/(windows-1):.2f}' #name of the end of the config
        file = 'top_eq_'+perc+'.pdb' #initial configurations
        directory = f'{int(float(perc)*100):02d}' #directories
        
        os.mkdir('bstates/'+directory)
        print('/'+directory)
        
        shutil.copyfile('config/'+file,'bstates/'+directory+'/top.pdb') #copy the file into the dir
        
        file = 'config/'+file

        ipos=get_initial_pos(file)

        prot=protein(ipos,args.n_chains)

        pos=CM(fasta_WT,prot)

        cl,centers=clust(pos,rc,L,2)

        print('Number of clusters: ',len(cl['frame 0']))
        rad = []
        for i in cl['frame 0']:
            print(f'cluster {i:} size: {cl["frame 0"][i]["size"]:} and radius: {cl["frame 0"][i]["rad"]:.6f}')
            
            rad.append(cl["frame 0"][i]["size"])   
        with open('bstates/'+directory+'/pcoord.init','w') as f:
            try:
                f.write(str(np.max(rad)))
                bins.append(np.max(rad))
            except ValueError:
                f.write(str(0))
                bins.append(0)
            pass
                
        with open('bstates/pcoord.init','a') as f:
            try:
                f.write(str(np.max(rad))+'\n')
            except ValueError:
                f.write(str(0)+'\n')
            pass

        with open('bstates/bstates.txt','a') as f:
            f.write(str(config)+' '+str(round(1/windows,2))+' '+directory+'\n')
                
        print('')
        
    #suggested bins
        
    bin_min = np.min(bins)
    bin_max = np.max(bins)+5.
    
    sug_bins=[]
    
    for i in range(windows):
        bin_lim = round((bin_max-bin_min)/(windows-1)*i,2)
        sug_bins.append(bin_lim)
    
    sug_bins.append('inf')
    
    string=''

    for i in sug_bins:
        
        if i == 'inf':
            string += i
        else:
            string += str(i)+','
    print('suggested bin distribution: ['+string+']') 
        
        
if __name__ == '__main__':
    directories(args.windows,args.seq,args.rc,args.L,args.n_chains)
    # for i in range(10):
    #     rc = rc_ini*(rc_fin-rc_ini)/10*i
    #     directories(args.windows,args.seq,args.rc,args.L)




