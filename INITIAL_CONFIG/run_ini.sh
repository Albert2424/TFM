#!/bin/bash

seq='WT' #protein name (WT or shuffle)
temp=300 #Temperature
cutoff=4.0 #Cutoff distance
wind=2 #number of windows desired for WESTPA (min 2)
n=100 #number of chains 
L=300 #cubic box size
steps=4000000 #steps for the simulation to reach equilibrium

read -p "This will remove all previous configurations, continue? (y/n) " answer

if [ $answer = 'y' ]; then

	if [ ! -d $seq ]; then
	    mkdir $seq  
	fi

	if [ ! -d $seq/$temp ]; then
	    mkdir $seq/$temp
	fi

	if [ ! -d 'config' ]; then
	    mkdir 'config'
	else
	    rm -r 'config'
	    mkdir 'config'
	fi

	python ./init_conf.py --windows $wind --n_chains $n --L $L --seq $seq
	python ./simulate_ini.py --seq $seq --temp $temp --cutoff $cutoff --windows $wind --n_chains $n --steps $steps

	rm -r $seq

else
    echo 'Stopped initial configuration generation'	
fi


