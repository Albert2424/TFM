#!/bin/bash
#SBATCH -J WESTPA
#SBATCH -e WEST.err
#SBATCH -o WEST.out
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=64
#SBATCH -N 1
##SBATCH --ntasks-per-node=25
#SBATCH --time=0-01:00:00
#SBATCH --mem=48G
#

module load cesga/2020 gcc/system openmpi/4.0.5_ft3 openmm/7.5.1-cuda-system mdanalysis

rm -f *.log
cd common_files
python analyse.py #generate pickle of the proteins.
cd ..
date

if [ ! -d "traj_segs" ]; then
    ./init.sh
fi

./run.sh
date
