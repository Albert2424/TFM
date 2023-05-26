#!/bin/bash

# Set up environment for westpa
# Actviate a conda environment containing westpa, openmm and mdtraj;
# you may need to create this first (see install instructions)

module load cesga/2020 gcc/system openmpi/4.0.5_ft3 openmm/7.5.1-cuda-system mdanalysis

#conda activate westpa2
export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)

export SWROOT=/mnt/netapp1/Optcesga_FT2_RHEL7/2020/gentoo/22072020

export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=300
export BASH=$SWROOT/bin/bash
export PERL=$SWROOT/usr/bin/perl
export ZSH=$SWROOT/bin/zsh
export IFCONFIG=$SWROOT/bin/ifconfig
export CUT=$SWROOT/usr/bin/cut
export TR=$SWROOT/usr/bin/tr
export LN=$SWROOT/bin/ln
export CP=$SWROOT/bin/cp
export RM=$SWROOT/bin/rm
export SED=$SWROOT/bin/sed
export CAT=$SWROOT/bin/cat
export HEAD=$SWROOT/bin/head
export TAR=$SWROOT/bin/tar
export AWK=$SWROOT/usr/bin/awk
export PASTE=$SWROOT/usr/bin/paste
export GREP=$SWROOT/bin/grep
export SORT=$SWROOT/usr/bin/sort
export UNIQ=$SWROOT/usr/bin/uniq
export HEAD=$SWROOT/usr/bin/head
export MKDIR=$SWROOT/bin/mkdir
export ECHO=$SWROOT/bin/echo
export DATE=$SWROOT/bin/date
