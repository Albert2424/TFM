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

def CM(fasta,prot,L):
    
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
            
            #PBC TO CALCULATE THE CENTER OF MASS OF THE CLUSTER
            X = np.array(p)[:,0]
            Y = np.array(p)[:,1] 
            Z = np.array(p)[:,2]
            
            theta_x = X/L*2*np.pi
            theta_y = Y/L*2*np.pi
            theta_z = Z/L*2*np.pi
            
            chi_x = np.cos(theta_x)
            zita_x = np.sin(theta_x)
            chi_y = np.cos(theta_y)
            zita_y = np.sin(theta_y)
            chi_z = np.cos(theta_z)
            zita_z = np.sin(theta_z)
            
            av_chi_x = np.sum(np.dot(chi_x,mass))/M
            av_zita_x = np.sum(np.dot(zita_x,mass))/M
            av_chi_y = np.sum(np.dot(chi_y,mass))/M
            av_zita_y = np.sum(np.dot(zita_y,mass))/M
            av_chi_z = np.sum(np.dot(chi_z,mass))/M
            av_zita_z = np.sum(np.dot(zita_z,mass))/M 
            
            av_theta_x = np.arctan2(-av_zita_x,-av_chi_x) + np.pi
            av_theta_y = np.arctan2(-av_zita_y,-av_chi_y) + np.pi
            av_theta_z = np.arctan2(-av_zita_z,-av_chi_z) + np.pi
            
            x_c = L*av_theta_x/(2*np.pi) 
            y_c = L*av_theta_y/(2*np.pi) 
            z_c = L*av_theta_z/(2*np.pi)
            
            r_CM[frame].append([x_c,y_c,z_c])
        
    return r_CM

def get_points(fasta,prot):

    """
    Gets 10 points of every chain for every frame.
    -------------------------------------------------------------
    
    INPUT:
    -------
        fasta: integer. Number of particles per chain.               
        
        prot: dictionary containing all frames of the trajectory.
   
    OUTPUT:
    --------
        r_prot: dictionary containing for all frames the ten positions of 
        each chain.                 
    
    """
    
    r_prot = {}
    mask = [i for i in range(0,len(fasta),len(fasta)//9)]
    for frame in prot:
        r_prot[frame] = []
        for p in prot[frame]:
            r_prot[frame].append(p[mask])
    return r_prot

def clust_detection(X, radius, L,min_samples=2):
    
    """
    Uses the DBSCAN to sample the different clusters of the system.
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
    
    dr = X[:,None, :, None,:]-X[None, :, None,:,:]
    dr = np.abs(dr - np.rint(dr/L)*L)#apply PBC

    dists = np.sqrt(np.sum((dr)**2, axis=-1))
    dists = dists.min(axis=(-1,-2))
    
    # print(np.min(dists[dists!=0])) #prints the minimum dist. rc around this value?
    
    db = DBSCAN(eps=radius,metric='precomputed', min_samples=min_samples).fit(dists)
    labels = db.labels_
    
    return labels

def get_radius(clusters,fasta,prot,frame,L):
    
    """
    Computes the radius of the cluster for a certain frame.
    -------------------------------------------------------------
    
    INPUT:
    -------
        clusters: dictionary containing every cluster of the frame. Each cluster 
                  has its own dictionary that contains the chains that are part of it, the
                  position of the center of mass of said chains and the size and radius
                  of the cluster:
            
            clusters = {'frame 0':{0:{'chains':[0,1,2...],
                                     'pos':[np.array(3),np.array(3),..],
                                     'size': 14},
                                   
                                   1:{...}
                                   
                                   .
                                   .
                                   .
                                   }
                        .
                        .
                        .
                        }
            
        fasta: integer. Number of particles per chain.               
        
        prot: dictionary containing all frames of the trajectory.
        
        frame: str. Name of the current frame: 'frame 0', 'frame 1'...
   
    OUTPUT:
    --------
        Adds to clusters the radius of every frame so now clusters is like:
            
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

    """
    
    pos=CM(fasta,prot,L) #get the center of mass of every chain
    
    # Get the centers of each cluster

    centers = []
    for clust in clusters[frame]:
        chain_per_clust = clusters[frame][clust]['chains'] #chains of the clust 
        positions = np.array(pos[frame])[chain_per_clust] #CM of the clust
        
        #PBC TO CALCULATE THE GEOMETRIC CENTER OF THE CLUSTER
        X = np.array(positions)[:,0]
        Y = np.array(positions)[:,1] 
        Z = np.array(positions)[:,2]
        
        theta_x = X/L*2*np.pi
        theta_y = Y/L*2*np.pi
        theta_z = Z/L*2*np.pi
        
        chi_x = np.cos(theta_x)
        zita_x = np.sin(theta_x)
        chi_y = np.cos(theta_y)
        zita_y = np.sin(theta_y)
        chi_z = np.cos(theta_z)
        zita_z = np.sin(theta_z)
        
        av_chi_x = np.mean(chi_x)
        av_zita_x = np.mean(zita_x) 
        av_chi_y = np.mean(chi_y)
        av_zita_y = np.mean(zita_y) 
        av_chi_z = np.mean(chi_z)
        av_zita_z = np.mean(zita_z) 
        
        av_theta_x = np.arctan2(-av_zita_x,-av_chi_x) + np.pi
        av_theta_y = np.arctan2(-av_zita_y,-av_chi_y) + np.pi
        av_theta_z = np.arctan2(-av_zita_z,-av_chi_z) + np.pi
        
        x_c = L*av_theta_x/(2*np.pi) 
        y_c = L*av_theta_y/(2*np.pi) 
        z_c = L*av_theta_z/(2*np.pi)
        
        centers.append([x_c,y_c,z_c])
          
        distances = np.linalg.norm(positions - centers[clust], axis=1)
        
        #rotation radius
        # clusters[frame][clust]['rad'] = np.max(distances)#add the radius 
        
        #hidrodynamic radius
        clusters[frame][clust]['rad'] = 1/(np.average(1/distances))
        clusters[frame][clust]['error'] = 1/(np.average(1/distances)**2)*np.std(1/distances)            


def clust(pos,dist,L,min_size,fasta,prot):
    
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
                     
        fasta: integer. Number of particles per chain.               
        
        prot: dictionary containing all frames of the trajectory.
   
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
        # get_radius(clusters,fasta,prot,frame,L) #adds radius to the dict
        
    return clusters


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
    
    print(f'Number of clusters in {frame:}: ',len(cl[frame]))
    rad = []
    for i in cl[frame]:
        print(f'cluster {i:} size: {cl[frame][i]["size"]:}')
        
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
    
    #clusters from the parent
    filename = "parent.pdb"
    ipos = get_initial_pos(filename)
    prot = protein(ipos,args.n_chains)
    pos = get_points(fasta_WT,prot)
    cl = clust(pos,args.rc,args.L,2,fasta_WT,prot)
    
    dist1 = generate_pcoord('frame 0',cl)
    
    #clusters from the trajectory
    filename = "seg.pdb"
    traject = 'seg.dcd'
    ipos = get_traj(traject,filename)
    prot = protein(ipos,args.n_chains)
    pos = get_points(fasta_WT,prot)
    cl = clust(pos,args.rc,args.L,2,fasta_WT,prot)
    
    
    dist2 = generate_pcoord('frame 0',cl)
    dist3 = generate_pcoord('frame 1',cl)
    
    
    d_arr = [dist1,dist2,dist3]       
    np.savetxt("dist.dat", d_arr) #save to dist.dat for westpa to get the pcoord
