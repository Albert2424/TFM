#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:31:16 2023

@author: albert
"""

from scipy.ndimage import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

col = cm.viridis(np.linspace(0, 1, 3))

residues = pd.read_csv('residues.csv')
residues.set_index('one', inplace=True)
seq='''RTEQLFSKNWSEQGGPNLERLNATNRDDSFSRYGTFQHRPSAGAPLHQMGETPAGQTDSTKSPT
PGESRIPLPGNVRIPGFPKKQNLSMIERHASNHQGNKNFHPFISSSFEVWEQESNWGADDATKQ
PGASGPGPPFYSPADHSWSMLCMQAQAGPISDHSVLPRYAFFHDHGLAIDVHAQLGQPNFATNN
DFQGGQFAQRSQGLPAAELSGNTHPESASNNSNPQPQNERGISSKPGGGYNSNNHHGWFPFAPK
TSLSAVIDQGLEYQSSPHSLSGAASPGALSQLKEDFSFSHFSIHPQRQNPIINPNHTGAQPANQ
HRDIGESIDGSGPGLESLHGPFSGSKRPHTQPAALNLAGVFIQENYATSHLSESGKLFRLHQKS
TSPPSHTASGSPHKKDPQEGFILGNTNAFIQQHRHGFLMNSSSPPLLWSPIHTRGSSDRPADPA'''.replace('\n', '').replace(' ', '')
hydrophob = [residues.loc[r].lambdas for r in seq] #seq és la sequencia en un string

plt.figure(figsize=(10,10)) 

plt.plot(gaussian_filter1d(hydrophob, 20), '-',label='shuffle', color=col[0])

seq="""
            MGDYG FGVLV QSNTG NKSAF PVRFH PHLQP PHHHQ NATPS PAAFI NNNTA ANGSS AGSAW
            LFPAP ATHNI QDEIL GSEKA KSQQQ EQQDP LEKQQ LSPSP GQEAG ILPET EKAKS EENQG
            DNSSE NGNGK EKIRI ESPVL TGFDY QEATG LGTST QPLTS SASSL TGFSN WSAAI APSSS
            TIINE DASFF HQGGV PAASA NNGAL LFQNF PHHVS PGFGG SFSPQ IGPLS QHHPH HPHFQ
            HHHSQ HQQQR RSPAS PHPPP FTHRN AAFNQ LPHLA NNLNK PPSPW SSYQS PSPTP SSSWS
            PGGGG YGGWG GSQGR DHRRG LNGGI TPLNS ISPLK KNFAS NHIQL QKYAR PSSAF APKSW
            MEDSL NRADN IFPFP DRPRT FDMHS LESSL IDIMR AENDT IKART YGRRR GQSSL FPMED
            GFLDD GRGDQ PLHSG LGSPH CFSHQ NGE
            """.replace('\n', '').replace(' ', '')
hydrophob = [residues.loc[r].lambdas for r in seq] #seq és la sequencia en un string


plt.plot(gaussian_filter1d(hydrophob, 20), '-',label='WT',color=col[1])
plt.xlabel('Residue',fontsize=20)
plt.ylabel('Hydrophobicity',fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.legend(fontsize=20)
plt.tight_layout()

plt.savefig('hydro.pdf')
plt.show()