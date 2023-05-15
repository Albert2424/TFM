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

parser = ArgumentParser()
parser.add_argument('--seq',nargs='?',const='', type=str)
parser.add_argument('--rc',nargs='?',const='', type=float)
parser.add_argument('--L',nargs='?',const='', type=float)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
args = parser.parse_args()
# name='WT'


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

def clust_detection(X, radius, L,min_samples=2):
    
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
        L: float. Size of the simulation box.
   
    OUTPUT:
    --------
        dbscan.labels_: 1D array containing the label of each particle. Particles
                        with the same label correspond to the same cluster. Particles
                        with -1 label are not part of any cluster.
    
    """
    X=np.array(X) #avoids [array(),array(),...]
    def my_pdist(x,y,L):
        dr = np.abs(x-y)
        dr = np.abs(dr - np.rint(dr/L)*L)
        return np.linalg.norm(dr)


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
   
    OUTPUT:
    --------
        clusters: dictionary containing every cluster of the frame. Each cluster 
                  has its own dictionary that contains the chains that are part of it, the
                  position of the center of mass of said chains and the size and radius
                  of the cluster:
            
            clusters = {'frame 0':{0:{'chains':[0,1,2...],
                                     'pos':[np.array(3),np.array(3),..],
                                     'size': 14,
                                     'radius': 1.453},
                                   
                                   1:{...}
                                   
                                   .
                                   .
                                   .
                                   }
                        .
                        .
                        .
                        }
            
        centers: array containing the xyz coordinates of every cluster center.
                        
    
    """
    
    clusters={}
    
    for frame in pos:
        # print(frame)
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
        
        for frame in clusters:
            centers = []
            for clust in clusters[frame]:
                centers.append([np.mean(np.array(clusters[frame][clust]['pos'])[:,0]),
                                np.mean(np.array(clusters[frame][clust]['pos'])[:,1]),
                                np.mean(np.array(clusters[frame][clust]['pos'])[:,2])])
                  
                distances = np.linalg.norm(np.array(clusters[frame][clust]['pos']) - centers[clust], axis=1)
                
                #rotation radius
                # clusters[frame][clust]['rad'] = np.max(distances)#add the radius 
                
                #hidrodynamic radius
                clusters[frame][clust]['rad'] = 1/(np.average(1/distances))
                    
    return clusters,np.array(centers)


def generate_pcoord(frame,cl):
    
    """
    Detects clusters formed in every frame of the simulation.
    -------------------------------------------------------------
    
    INPUT:
    -------
        frame: name of the frame that you want the pcoord to be obtained from.
               It should be a string like "frame 0", "frame 1", etc
        
        cl: dictionary containing every cluster of the frame. Each cluster 
            has its own dictionary that contains the chains that are part of it, the
            position of the center of mass of said chains and the size and radius
            of the cluster:
                
                clusters = {'frame 0':{0:{'chains':[0,1,2...],
                                         'pos':[np.array(3),np.array(3),..],
                                         'size': 14,
                                         'radius': 1.453},
                                       
                                       1:{...}
                                       
                                       .
                                       .
                                       .
                                       }
                            .
                            .
                            .
                            }
               
    OUTPUT:
        dist: pcoord obtained, corresponds to the radius of the bigger cluster.
    --------

    
    """
    
    print(f'Number of clusters in frame {frame:}: ',len(cl[frame]))
    rad = []
    for i in cl[frame]:
        print(f'cluster {i:} size: {cl[frame][i]["size"]:} and radius: {cl[frame][i]["rad"]:.6f}')
        
        rad.append(cl[frame][i]["size"])
    
    print('')
    try:
        dist = np.max(rad)
    except ValueError:
        dist = 0
        
    return dist


if __name__ == "__main__":
    proteins = initProteins()
    fasta_WT = proteins.loc[args.seq].fasta
    
    filename = "top.pdb"
    traject = 'traj.dcd'
    
    ipos=get_traj(traject,filename)
    
    prot=protein(ipos,args.n_chains)
    
    pos=CM(fasta_WT,prot)
    
    cl,centers=clust(pos,args.rc,args.L,2)
    
    dist1 = generate_pcoord('frame 0',cl)
    dist2 = generate_pcoord('frame 1',cl)
    dist3 = generate_pcoord('frame 2',cl)
    d_arr = [dist1,dist2,dist3]       
    np.savetxt("dist.dat", d_arr)


