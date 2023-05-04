#!/bin/bash

#read the simulation variables

seq='WT'       #name of the studied protein (str)
temp=320       #temperature (int)
cutoff=4.0     #cutoff (float)
tau=4000        #steps of the simulation (tau in westpa) (int)
n_chains=100   #number of chains of the system (int)

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/bstate.pdb .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/simulate.py > simulate.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/cluster.py > cluster.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/top.pdb parent.pdb
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/simulate.py > simulate.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/cluster.py > cluster.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/top.pdb parent.pdb
fi

# Run the dynamics with OpenMM
python simulate.py --name $seq --temp $temp --cutoff $cutoff --steps $tau --n_chains $n_chains

#Calculate pcoord with MDAnalysis
python analyse.py #initialize the functions on analyse to be used in cluster.py
python cluster.py --seq $seq > clust.log
cat dist.dat > $WEST_PCOORD_RETURN

# Clean up
rm -f *.py *.csv *.pkl dist.dat *.npy *.log 
rm -f parent.pdb bstate.pdb WT.dcd
rm -r __pycache__
