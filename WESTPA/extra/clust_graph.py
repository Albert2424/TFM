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

    """
    config = []
    with open('bstates_'+seq+'/bstates.txt','r') as f:
        for line in f:
            config.append(line.split()[-1])
    
    
    # for c in range(windows):
    #     perc = f'{int(round((c/(windows-1)*n_chains),0)):02d}' #name of the end of the config
    #     file = perc #initial configurations
    #     config.append(file)
    return config

def get_radius(positions,clust,frame,clusters):
    
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
    
        
          
    distances = np.linalg.norm(positions[0], axis=1)
    
    #rotation radius
    # clusters[frame][clust]['rad'] = np.max(distances)#add the radius 
    
    #hidrodynamic radius
    clusters[frame][clust]['rad'] = 1/(np.average(1/distances)) # 1/<r^(-1)>
    clusters[frame][clust]['error'] = 1/(np.average(1/distances)**2)\
        *np.sqrt(np.std(1/distances))/len(distances) # 1/<r^(-1)>^2 * sqrt(std(<r^(-1)>))/N

def get_clusts(clusters,fasta,prot,frame,L,t):
    
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
    
    # Get the centers of each cluster
    biggest = 0
    for clust in clusters[frame]:
        size = clusters[frame][clust]['size']

        if size > biggest:
            max_clust = clust
            biggest = size
            
    c = max_clust        
    chains = clusters[frame][c]['chains'] #chains of the clust 
    # positions = np.array(pos[frame])[chain_per_clust] #atoms of the clust
    
    #PBC TO CALCULATE THE GEOMETRIC CENTER OF THE CLUSTER
    cluster_resi = t.top.select("chainid "+' '.join([str(i) for i in chains]))
    t = t.atom_slice(cluster_resi)
    bonds=np.array([(i,i+1) for i in np.arange(t.n_atoms-1)], dtype=np.int32)
    t.make_molecules_whole(inplace=True,sorted_bonds=bonds)
    t.center_coordinates()
    positions = t.xyz
    
    center = [0,0,0]
    
    get_radius(positions,c,frame,clusters)
    
    rad = clusters[frame][c]['rad']
    e = clusters[frame][c]['error']
        
    return center,biggest,chains,rad,positions,e

def rel_dist(config,fasta,n_chains,L,seq):
    dist = {}
    rad = []
    err = []
    for c in config:
        print(c)
        filename = 'bstates_'+seq+'/'+c+'/top.pdb'
        ipos = get_initial_pos(filename)
        prot = protein(ipos,n_chains)
        pos = get_points(fasta,prot)
        cl = clust(pos,17.,L,2,fasta,prot)#get the clusters
        t = md.load(filename)
        
        center, size, chains, radius,pos,e = get_clusts(cl,fasta,prot,'frame 0',L,t)#compute the center
                                                                  #of the clust, its
                                                                  #size and chains that
        rad.append(radius)                                        #are part of it.
        err.append(e)
        dist[c] = []
        N = len(fasta)
        count=1
        for i in chains:
            i += 1
            dists = pos[0][(count-1)*N:(count)*N] - center #select chain
            # dists = np.abs(dr - np.rint(dr/L)*L)#apply PBC
            distances = np.linalg.norm(dists, axis=1)
            dist[c].append(distances) 
            count+=1

    return dist,rad,err


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
    
    
    data_seq = []
    plt.figure(figsize=(10,10)) 
    out = ['00']
    col = iter(cm.viridis(np.linspace(0, 1, len(dist)-len(out))))
    
    for i in dist:
        if i not in out:
            av = np.mean(dist[i],axis=0)
            plt.plot(range(len(fasta)),av,
                      label=f'n = {int(i):}',color=next(col),marker='o',
                      markersize=5,markeredgewidth=0.5,markeredgecolor='black')
            data_seq.append(av)
    # plt.xticks(range(len(fasta)), fasta)    
    plt.xlabel(seq,fontsize=20)
    plt.ylabel('r (nm)',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()

    plt.savefig('dist_graph_'+seq+'.pdf')
    plt.show()
    
    return data_seq
    

def plot_dens(dist,fasta,n_chains,rad,seq,err):
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
    
    
    
    plt.figure(figsize=(10,8)) 
    out = ['00','11']
    col = cm.viridis(np.linspace(0, 1, len(dist)-len(out)))

    count = 0
    for i in dist:
        dist_tot = np.asarray([])
        if i not in out:
            for j in dist[i]:
                dist_tot = np.concatenate((dist_tot, j))
            
            radius = rad[count]
            e = err[count]
            # e=1
            # print(e)
            # bins = [35/49*i for i in range(50)]
            bins= np.linspace(0, 35**3, 50)**(1/3)
            hist,bins=np.histogram(dist_tot,bins=bins)
                
            hist = np.array(hist,dtype='float64')
            hist /=  4/3*np.pi*(bins[1:]**3-bins[:-1]**3)
            rad_c = (bins[1:]+bins[:-1])/2

            plt.plot(rad_c,hist, color=col[count-len(out)],label='n = '+ str(int(int(i)/100*n_chains)))
            plt.vlines(radius,0,2.5,colors=col[count-len(out)],linestyles='dashed')
            plt.axvspan(radius-e, radius+e, alpha=0.3, color=col[count-len(out)])
        count += 1
        
    # plt.xticks(range(len(fasta)), fasta)  
    plt.title(seq,fontsize=25)
    plt.xlabel('r (nm)',fontsize=20)
    plt.ylabel(r'$\rho \: (residues/nm^3)$',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.ylim((0,2.5))
    plt.xlim((4.5,35))
    plt.tight_layout()
    plt.savefig('dens_graph_'+seq+'.pdf')
    plt.show()
    
    
    n = [int(i) for i in dist]
    
    plt.figure()
    plt.scatter(n[1:],rad[1:])
    plt.plot(n[1:],3*np.log(n[1:]))
    plt.show()
        
            
#%%     
if __name__ == '__main__':
    seq = 'WT'
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    n_chains = 250
    L = 339
    config = config_list(10,n_chains,seq)
    # config=['00','11','22','33','44','56','67','78','89','100']
    
    print(f'analysing {seq}...\n')
    dist,rad,err = rel_dist(config,fasta,n_chains,L,seq)
    data_WT = plot(dist, fasta,n_chains,seq)
    plot_dens(dist, fasta,n_chains,rad,seq,err)

    seq = 'shuffle'
    proteins = initProteins()
    fasta = proteins.loc[seq].fasta
    
    print(f'analysing {seq}...\n')
    dist,rad,err = rel_dist(config,fasta,n_chains,L,seq)
    data_shuffle = plot(dist, fasta,n_chains,seq)
    plot_dens(dist, fasta,n_chains,rad,seq,err)

#%%  
    #COMPARISON PLOT
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    N = len(fasta)
    col1 = iter(cm.viridis(np.linspace(0, 1, 8)))
    col2 = iter(cm.viridis(np.linspace(0, 1, 8)))
    
    config = config_list(10,n_chians,seq)
    config = [int(i) for i in config[2:]]
    config = np.array(np.array(config)/100*n_chains,dtype='int')
    ticks = [(i)/(np.max(config)) for i in config]
    ticks = [i-ticks[0]*(1-ticks.index(i)/len(ticks)) for i in ticks]
    
    fig, axs = plt.subplots(1,2,figsize=(12,6),sharey=True)
    
    #shuffle
    plt.subplot(121)
    plt.title('Shuffle',fontsize=25)
    for i in data_shuffle:
        plt.plot(range(N),i,
                  color=next(col1),marker='o',
                  markersize=5,markeredgewidth=0.000,markeredgecolor='black')
    
    
    plt.xlabel('Residues',fontsize=20)
    plt.ylabel('r (nm)',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    
    #WT
    plt.subplot(122)
    plt.title('WT',fontsize=25)
    for i in data_WT:
        plt.plot(range(N),i,
                  color=next(col2),marker='o',
                  markersize=5,markeredgewidth=0.000,markeredgecolor='black')
    
    plt.xlabel('Residues',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    
    #colorbar

    divider = make_axes_locatable(axs[1])
    ax_cb = divider.append_axes("right", size="7%", pad="5%")
    cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.viridis, orientation='vertical')
    cb1.ax.get_yaxis().set_ticks(ticks)
    cb1.ax.set_yticklabels(config)
    cb1.set_label('Chains in the cluster', rotation=270, fontsize=20, labelpad=20)
    plt.gcf().add_axes(ax_cb)
    
    
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.tight_layout()
    
    plt.savefig('comparison.pdf')
    plt.show()
    
#%%

a = np.linspace(5, 70, 33)
print(list(a))
