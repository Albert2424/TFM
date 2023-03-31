#!/bin/bash

#read the simulation variables

name='WT'      #name of the studied protein (str)
temp=320       #temperature (int)
cutoff=4.0     #cutoff (float)
tau=400        #steps of the simulation (tau in westpa) (int)
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
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/top.pdb ./parent.pdb
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/simulate.py > simulate.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/top.pdb ./parent.pdb
fi

# Run the dynamics with OpenMM
python simulate.py --name $name --temp $temp --cutoff $cutoff --steps $tau --n_chains $n_chains

#Calculate pcoord with MDAnalysis
python $WEST_SIM_ROOT/common_files/cluster.py --name $name > clust.log
cat dist.dat > $WEST_PCOORD_RETURN

cp top.pdb $WEST_RESTART_RETURN/parent.pdb
# Clean up
rm -f *.py *.csv *.pkl dist.dat *.npy
rm -r __pycache__
