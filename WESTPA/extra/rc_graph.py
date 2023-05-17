#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 13:44:53 2023

@author: albert

rc graph
"""
from cluster import *
from analyse import *
from argparse import ArgumentParser
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN,KMeans
import matplotlib.cm as cm


parser = ArgumentParser()
parser.add_argument('--seq',nargs='?',const='', type=str)
parser.add_argument('--windows',nargs='?',const='', type=int)
parser.add_argument('--L',nargs='?',const='', type=float)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
args = parser.parse_args()



def config_list(windows):
    config = []
    expected = []
    for c in range(windows):
        perc = f'{c/(windows-1):.2f}' #name of the end of the config
        file = 'config/top_eq_'+perc+'.pdb' #initial configurations
        config.append(file)
        expected.append(int(c/(windows-1)*100))
    return config,expected
       
def rc_graph(config,proteins,fasta,n_chains,L):
    clust_list={}
    for i in range(max_it):
        rc = rc_ini+(rc_fin-rc_ini)/(max_it-1)*i
        print(f'~~~~ Current rc: {rc:.2f} ({i+1:}/{max_it:})~~~~')
        print('')
        clust_list[rc] = []
        for c in config:
            print(c)
            filename = c
            ipos=get_initial_pos(filename)
            prot=protein(ipos,n_chains)
            pos=CM(fasta_WT,prot)
            cl,centers=clust(pos,rc,L,2)
            aux = []
            for i in cl["frame 0"]:
                aux.append(cl["frame 0"][i]["size"])
            
            try:
                clust_list[rc].append(np.max(aux))
            except ValueError:
                clust_list[rc].append(0)
                
        with open('clust_rc.dat','a') as file:
            for i in clust_list[rc]:
                file.write(str(i)+'\n')
            file.write('\n')
        print('')
    return clust_list
        
            
def plot(clust_list,expected,max_it):
    print('plotting RC graph...')
    plt.figure(figsize=(10,10)) 
    col = iter(cm.viridis(np.linspace(0, 1, max_it)))
    for i in clust_list:
        # plt.plot(expected,np.array(clust_list[i])-np.array(expected),
        #          label=f'rc = {i:.2f}',color=next(col),marker='o',
        #          markersize=5, markeredgewidth=2.)
        plt.plot(expected,np.array(clust_list[i])-np.array(expected),
                 label='rc = '+i,color=next(col),marker='o',
                 markersize=5, markeredgewidth=0.5, markeredgecolor='black')
    
    plt.xlabel('ECS',fontsize=20)
    plt.ylabel('ECS-CS',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()

    plt.savefig('rc_graph.pdf')
    plt.show()
    
    
def read(filename,rc_ini,rc_fin,max_it):
    rc_ = {f'{rc_ini:.2f}':[]}
    with open(filename,'r') as file:
        count = 0
        rc = rc_ini+(rc_fin-rc_ini)/(max_it-1)*count
        for line in file:
            s = line.split()
            if len(s)>0:
                rc_[f'{rc:.2f}'].append(int(s[0])) 
            else:
                count += 1
                rc = rc_ini+(rc_fin-rc_ini)/(max_it-1)*count
                rc_[f'{rc:.2f}']=[]
    rc_.popitem()
    return rc_
                
                
            
    
if __name__ == '__main__':
    rc_ini = 20
    rc_fin = 45
    max_it = 10
    windows = 15
    
    config, expected = config_list(windows)
    proteins = initProteins()
    fasta_WT = proteins.loc[args.seq].fasta
    
    clust_list = rc_graph(config, proteins, fasta_WT, args.n_chains, args.L)
    # clust_list = read('clust_rc.dat',rc_ini,rc_fin,max_it)
    plot(clust_list,expected,max_it)
