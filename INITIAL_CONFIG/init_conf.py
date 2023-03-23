#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:50:47 2023

@author: albert
"""

"""INITIAL CONDITIONS FOR WESTPA"""

import simtk.openmm as openmm
from simtk.openmm import unit
from simtk.openmm import app
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import random
from analyse_ini import *
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--windows',nargs='?',const='', type=int)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
parser.add_argument('--L',nargs='?',const='', type=int)
parser.add_argument('--seq', nargs='?',const='', type=str)
args = parser.parse_args()

# windows=3
# n_chains=100
# L=300
# seq='WT'


prot={'WT':"""
            MGDYG FGVLV QSNTG NKSAF PVRFH PHLQP PHHHQ NATPS PAAFI NNNTA ANGSS AGSAW
            LFPAP ATHNI QDEIL GSEKA KSQQQ EQQDP LEKQQ LSPSP GQEAG ILPET EKAKS EENQG
            DNSSE NGNGK EKIRI ESPVL TGFDY QEATG LGTST QPLTS SASSL TGFSN WSAAI APSSS
            TIINE DASFF HQGGV PAASA NNGAL LFQNF PHHVS PGFGG SFSPQ IGPLS QHHPH HPHFQ
            HHHSQ HQQQR RSPAS PHPPP FTHRN AAFNQ LPHLA NNLNK PPSPW SSYQS PSPTP SSSWS
            PGGGG YGGWG GSQGR DHRRG LNGGI TPLNS ISPLK KNFAS NHIQL QKYAR PSSAF APKSW
            MEDSL NRADN IFPFP DRPRT FDMHS LESSL IDIMR AENDT IKART YGRRR GQSSL FPMED
            GFLDD GRGDQ PLHSG LGSPH CFSHQ NGE
            """.replace('\n', '').replace(' ', ''),
            
     'Shuffle': """
                MGDYGFGVLVQSNTGNKSAFPVRFHPHLQPPHHHQNATPSPAAFINNNTAANGSSAGSAWLFP
                APATHNIQDEILGSEKAKSQQQEQQDPLEKQQLSPSPGQEAGILPETEKAKSEENQGDNSSEN
                GNGKEKIRIESPVLTGFDYQEATGLGTSTQPLTSSASSLTGFSNWSAAIAPSSSTIINEDASF
                FHQGGVPAASANNGALLFQNFPHHVSPGFGGSFSPQIGPLSQHHPHHPHFQHHHSQHQQQRRS
                PASPHPPPFTHRNAAFNQLPHLANNLNKPPSPWSSYQSPSPTPSSSWSPGGGGYGGWGGSQGR
                DHRRGLNGGITPLNSISPLKKNFASNHIQLQKYARPSSAFAPKSWMEDSLNRADNIFPFPDRP
                RTFDMHSLESSLIDIMRAENDTIKARTYGRRRGQSSLFPMEDGFLDDGRGDQPLHSGLGSPHC
                FSHQNGE
                """.replace('\n', '').replace(' ', '')     
     }



def initial_config(L,margin,n_chains,perc,N):
    
    """
    Provides the position of the center of the chains in the box.
    It takes into acount how many chains must be close together to
    induce cluster formation.
    -------------------------------------------------------------
    
    INPUT:
    -------
        L: float type. Length of the cubic box.
        
        margin: float type. Margin of the box where chains won't be placed 
                any further from.
                
        n_chains: integer type. Number of chains to be placed.
        
        perc: float type. (Between 0 and 1) indicates the percentge of chains 
             placed on the center of the box.
             
        N: integer type. Lenght of the chains.
        
    OUTPUT:
    --------
        xyz_clust: (int(n_chain*perc),3) float array type. Array with the xyz coordinates 
                   of the particles placed to induce the cluster.
                   
        xyz_non_clust: (int(n_chain*(1-perc)),3) float array type. Array with the xyz coordinates 
                   of the particles placed randomly along the box.                   

    """
    
    xyz = np.empty((n_chains,3)) #generates empty coordinate matrix
    
    L = L-margin
    
    #CLUSTER CHAINS
    
    #random xy distribution
    n_clust = int(n_chains*perc) #chains in cluster
    
    L_clust = int(np.sqrt(n_clust)*2) #limited space to induce a cluster.
    
    xyz_clust = np.empty((n_clust,3))
    
    if n_clust>0:
            
        # Create a list of all possible coordinate pairs within the lattice
        xy_coord = [[-L_clust/2+x, -L_clust/2+y] for x in range(L_clust+1) for y in range(L_clust+1)]
        # Shuffle the list
        random.shuffle(xy_coord)
        
        # Assign each number in sequence to a coordinate pair
        for i in range(n_clust):
            xyz_clust[i,0:2] = xy_coord[i]
            xyz_clust[i,-1] = 0.        
    
    #NON-CLUSTER CHAINS
    
    n = n_chains-n_clust
    xyz_non_clust = np.empty((n,3))
    xy_coord = []
        
    for i in range(L+1):
        for j in range(L+1):
                
            x = -L/2 + i 
            y = -L/2 + j
            
            if (x > L_clust or x < -L_clust) and (y > L_clust or y < -L_clust):
                xy_coord.append([x,y])
                
    # Shuffle the list
    random.shuffle(xy_coord)
    
    # Assign each number in sequence to a coordinate pair
    for i in range(n):
        xyz_non_clust[i,0:2] = xy_coord[i]
        xyz_non_clust[i,-1] = random.uniform(int(-L/2), int(L/2))
    
    return xyz_clust,xyz_non_clust

def topo(xyz,perc,n_chains,L,N):
    
    """
    Adds topology to the center of the different chains in order to create
    a .pdb file with the initial conditions.
    ----------------------------------------------------------------------
    
    INPUT:
    -------
        xyz: (n,3) float array type. n is the number of chains of the system. 
             Contains the xyz coordinates of the centers of the chains.
        
        per: float type. (Between 0 and 1) indicates the percentge of chains 
             placed on the center of the box.
             
        n_chains: integer type. Number of chains of the system.
        
        L: float type. Length of the cubic box.
        
        N: integer type. Lenght of the chains.
             
    OUTPUT:
    --------
        Returns a .pdb file named top_{per}.pdb with the initial configuration. 
            
    """
    
    import mdtraj as md
    
    pdb_file = 'config/top_{:4.2f}.pdb'.format(perc)
    
    top = md.Topology()
    pos = []
    for x,y,z in xyz:
        chain = top.add_chain()
        pos.append([[x,y,z+(i-N/2.)*.38] for i in range(N)])
        for resname in fasta:
            residue = top.add_residue(resname, chain)
            top.add_atom(resname, element=md.element.carbon, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))
            # if i==int(N/2):
                # print(top.atom(i))
            
    md.Trajectory(np.array(pos).reshape(n_chains*N,3), top, 0, [L,L,L], [90,90,90]).save_pdb(pdb_file)
    pdb = app.pdbfile.PDBFile(pdb_file)
    
    #change beta factor of the chains for visual representation.
    b=[]
    for i in range(n_chains):
        for j in range(N):
            b.append(j/N)
            
    traj=md.load_pdb(pdb_file)
    traj.save_pdb(filename=pdb_file,bfactors=b)
    
def west_conf(windows,L,margin,n_chains,N):
    
    """
    Generates a initial configuration with an induced cluster
    of different size for each window.
    -------------------------------------------------------------
    
    INPUT:
    -------
        windows: integer type. Number of windows.
    
        L: float type. Length of the cubic box.
        
        margin: float type. Margin of the box where chains won't be placed 
                any further from.
                
        n_chains: integer type. Number of chains to be placed.
             
        N: integer type. Lenght of the chains.
        
    OUTPUT:
    --------
        Files with the initial configurations for each window of the type 
        top_{per}.pdb                   

    """
    
    for i in range(windows):
        
        perc=i/(windows-1)
        
        # print(f'Generating config {perc:}...')
        
        clust, non_clust=initial_config(L,margin,n_chains,perc,N)
        # print('clust:',clust,'\n\n','non clust',non_clust)
        xyz=np.vstack((clust,non_clust)) #stack both configurations
        topo(xyz,perc,n_chains,L,N)
    
     
    print(f'{i+1:} configurations generated')
        
    

#%%

print('Generating configurations...')
proteins = initProteins()
fasta = proteins.loc[args.seq].fasta
# fasta = prot[args.seq]
  
N = len(fasta)
margin = 8

west_conf(args.windows,args.L,margin,args.n_chains,N)

# %%

# print('Generating configurations...')

# proteins = initProteins()
# fasta = proteins.loc[seq].fasta
# perc=0.3
  
# N = len(fasta)
# margin = 8

# west_conf(8,300,margin,100,N)


