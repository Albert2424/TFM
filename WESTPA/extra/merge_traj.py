#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:29:46 2023

@author: albert
"""

import os
import shutil
import mdtraj as md
import re
import numpy as np

'''
REQUIRES:
    - seg_logs of the <iteration> (ex: segs_logs)
    - seg_traj/<iteration> (ex: 000081)
'''

#Sets the order in which the files must be read.

sizes = [] 
for root, dirs, files in os.walk('/home/albert/TFM/RESULTS_ft3/seg_logs'):

    for file in files: 
        aux = []
        with open('seg_logs/'+file, "r") as f:
            for line in f:
                s = line.split()
                if line.startswith('cluster'):
                    aux.append(int(s[3]))
        sizes.append(f'{np.max(aux):03d}')
                    
    arr = np.vstack((files,sizes))   
    arr = arr[:,arr[1,:].argsort()] 
#generate the merged dcd
segs = len(files)

count = 1
for i in arr[0]:
    directory = i.split('-')[0]+'/'+i.split('-')[-1][:6]
    # print(directory)
    print(count,f' / {segs:}')
    
    if list(arr[0]).index(i) == 0:
        merged_traj = md.load_dcd(directory+'/seg.dcd', top=directory+'/seg.pdb')
    else:
        traj = md.load_dcd(directory+'/seg.dcd', top=directory+'/seg.pdb')
        merged_traj = merged_traj.join(traj)
    count += 1
        
merged_traj.save_dcd("merged_trajectory.dcd")

         
                  