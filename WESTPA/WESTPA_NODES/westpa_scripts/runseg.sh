#!/bin/bash

#read the simulation variables

seq='WT'       #name of the studied protein (str)
temp=320       #temperature (int)
cutoff=4.0     #cutoff (float)
tau=600000     #steps of the simulation (tau in westpa) (int)
n_chains=100   #number of chains of the system (int)
L=300.         #length of the box
rc=30.         #minimum radius to consider a chain is inside a cluster

#if [ -n "$SEG_DEBUG" ] ; then
#  set -x
#  env | sort
#fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

#ln -sv $WEST_SIM_ROOT/common_files/bstate.pdb .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/simulate.py > simulate.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/cluster.py > cluster.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/seg.pdb ./parent.pdb
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/simulate.py > simulate.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/cluster.py > cluster.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/analyse.py > analyse.py
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/residues.csv > residues.csv
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/proteins.pkl > proteins.pkl
  ln -sv $WEST_PARENT_DATA_REF/top.pdb ./parent.pdb
fi

export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES
echo "node name" $SLURMD_NODENAME

# Run the dynamics with OpenMM
python simulate.py --name $seq --temp $temp --cutoff $cutoff --steps $tau --n_chains $n_chains --L $L

#Calculate pcoord with MDAnalysis
python analyse.py #initialize the functions on analyse to be used in cluster.py
python cluster.py --seq $seq --rc $rc --L $L --n_chains $n_chains > clust.log
cat dist.dat > $WEST_PCOORD_RETURN

# Clean up
rm -f *.py *.csv *.pkl dist.dat *.npy *.log 
rm -f parent.pdb
rm -r __pycache__
