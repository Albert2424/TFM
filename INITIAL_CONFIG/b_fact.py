#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 12:34:15 2023

@author: albert

Change beta factor in the pbd files
"""

import mdtraj as md
n_chains=100
fasta=448
file= "config/top_eq_0.60.pdb"


b=[]
for i in range(n_chains):
    for j in range(fasta):
        b.append(j/fasta)
        
traj=md.load_pdb(file)
traj.save_pdb(filename=file,bfactors=b)


