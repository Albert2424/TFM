#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:33:14 2023

@author: albert

cluster homogeniety
"""

from cluster import *
from analyse import *
from argparse import ArgumentParser
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN,KMeans
import matplotlib.cm as cm
import seaborn as sns


parser = ArgumentParser()
parser.add_argument('--seq',nargs='?',const='', type=str)
parser.add_argument('--windows',nargs='?',const='', type=int)
parser.add_argument('--L',nargs='?',const='', type=float)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
args = parser.parse_args()

def config_list(windows):
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
        perc = f'{int(round((100*c/(windows-1)),0)):02d}' #name of the end of the config
        file = perc #initial configurations
        config.append(file)
    return config

def get_clusts(clusters,fasta,prot,frame,L):
    
    """
    Computes the center of every cluster, the size of the biggest
    cluster and the chains that are part of it.
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
        center: array (3). xyz coords of the center of the biggest cluster
        
        biggest: int. size of the biggest cluster
        
        chain_per_clust: list. chains of the biggest cluster.

    """
    
    pos=CM(fasta,prot,L) #get the center of mass of every chain
    
    # Get the centers of each cluster
    biggest = 0
    for clust in clusters[frame]:
        size = clusters[frame][clust]['size']

        if size > biggest:
            chain_per_clust = clusters[frame][clust]['chains'] #chains of the clust 
            positions = np.array(pos[frame])[chain_per_clust] #atoms of the clust
            
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

            center = [x_c,y_c,z_c]
            # print(center)
            biggest = size
            
            rad = clusters[frame][clust]['rad']
        
    return center,biggest,chain_per_clust,rad

def rel_dist(config,fasta,n_chains,L):
    dist = {}
    rad = []
    for c in config:
        print(c)
        filename = 'bstates/'+c+'/top.pdb'
        ipos = get_initial_pos(filename)
        prot = protein(ipos,n_chains)
        pos = get_points(fasta,prot)
        cl = clust(pos,17.,L,2,fasta,prot)#get the clusters
        
        center, size, chains, radius = get_clusts(cl,fasta,prot,'frame 0',L)#compute the center
                                                                  #of the clust, its
                                                                  #size and chains that
        rad.append(radius)                                                          #are part of it.
        dist[c] = []
        N = len(fasta)
        for i in chains:
            i += 1
            dr = ipos[0][(i-1)*N:(i)*N] - center #select chain
            dists = np.abs(dr - np.rint(dr/L)*L)#apply PBC
            distances = np.linalg.norm(dists, axis=1)
            dist[c].append(distances) 
            # print(len(chains),len(distances))
    return dist,rad


def plot(dist,fasta,n_chains,seq):
    """
    Plots the graph of the average distance of every residue to the 
    center of the cluster.
    -------------------------------------------------------------
    
    INPUT:
    -------
        dist: dict. contains the distances to the center for every chain of the cluster.
                    the keys are the different configurations while the
                    values are arrays containing the values for every configuration.
        
        fasta: list. fasta of the protein.
        
        n_chains: int. chains of the system.
        
    OUTPUT:
    --------
        Generates a plot that is saved as 'dist_graph.pdf'.

    """
    
    print('plotting distances graph...')
    
    
    
    plt.figure(figsize=(10,10)) 
    out = ['00','11']
    col = iter(cm.viridis(np.linspace(0, 1, len(dist)-len(out))))
    
    for i in dist:
        if i not in out:
            av = np.mean(dist[i],axis=0)
            plt.plot(range(len(fasta)),av,
                      label=f'n = {int(int(i)/100*n_chains):}',color=next(col),marker='o',
                      markersize=5,markeredgewidth=0.5,markeredgecolor='black')
        
    # plt.xticks(range(len(fasta)), fasta)    
    plt.xlabel(seq,fontsize=20)
    plt.ylabel('r (nm)',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()

    plt.savefig('dist_graph_'+seq+'.pdf')
    plt.show()
    

def plot_dens(dist,fasta,n_chains,rad,seq):
    """
    Plots the graph of the average distance of every residue to the 
    center of the cluster.
    -------------------------------------------------------------
    
    INPUT:
    -------
        dist: dict. contains the distances to the center for every chain of the cluster.
                    the keys are the different configurations while the
                    values are arrays containing the values for every configuration.
        
        fasta: list. fasta of the protein.
        
        n_chains: int. chains of the system.
        
    OUTPUT:
    --------
        Generates a plot that is saved as 'dist_graph.pdf'.

    """
    
    print('plotting density graph...')
    
    
    
    plt.figure(figsize=(10,10)) 
    out = ['00','11']
    col = cm.viridis(np.linspace(0, 1, len(dist)-len(out)))

    dist_tot = np.asarray([])
    count = 0
    for i in dist:
        if i not in out:
            for j in dist[i]:
                dist_tot = np.concatenate((dist_tot, j))
            
            radius = rad[count]
            hist,bins=np.histogram(dist_tot,bins=50)
            # dr = bins[1]-bins[0]
            # for k in range(len(hist)):
            #     # dv = 4/3*np.pi*dr**3*(3*count**2+3*count+1) #volum of a layer
            #     hist[k] /= dv
                
            hist = np.array(hist,dtype='float64')
            hist /=  4/3*np.pi*(bins[1:]**3-bins[:-1]**3)
            rad_c = (bins[1:]+bins[:-1])/2
            
            plt.plot(rad_c,hist, color=col[count-len(out)],label=f'n = {int(int(i)/100*n_chains):}')
            plt.vlines(radius,0,16,colors=col[count-len(out)],linestyles='dashed')
        count += 1
        
    # plt.xticks(range(len(fasta)), fasta)    
    plt.xlabel(seq,fontsize=20)
    plt.ylabel(r'$\rho$ (nm)',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()
    # plt.ylim((0,0.009))
    # plt.xlim((0,40))

    plt.savefig('dens_graph_'+seq+'.pdf')
    plt.show()
        
            
#%%     
if __name__ == '__main__':
    config = config_list(10)
    seq = 'WT'
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    n_chains = 150
    L = 343
    
    dist,rad = rel_dist(config,fasta,n_chains,L)
    plot(dist, fasta,n_chains,seq)
    plot_dens(dist, fasta,n_chains,rad,seq)
  
