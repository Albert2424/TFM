# INITIAL CONFIGURATION

In order to start the dynamics from different starting conditions, some specific initial distributions are generated. The aim of generating distributions with different aggregate sizes is to optimize the WESTPA simulation of the dynamics of aggregation. 

To generate the initial configurations we need three codes: `init_config.py`, which provides a first spatial distribution with stiff protein distribution and `simulate_ini.py` and `analyse_ini.py` which are meant to bring the configuration to equilibrium while maintaining the chains that are far away from the cluster fixed. Both `simulate_ini.py` and `analyse_ini.py` are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper _Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range_. See the original code in **[here]**.



[here]: https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code


## REQUIREMENTS

The code has been created using _Python 3.10.9_. Other important packages needed are:

* OpenMM ([8.0.0])
* MDAnalysis ([2.4.2]) 
* MDtraj ([1.9.7])

[8.0.0]: http://docs.openmm.org/7.0.0/userguide/application.html
[2.4.2]: https://www.mdanalysis.org/pages/installation_quick_start/
[1.9.7]: https://www.mdtraj.org/1.9.7/installation.html

All packages may be installed using [Miniconda] with the `pip` command of the conda prompt.

[Miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html


## EXECUTION

In order to execute the program you can both use a local machine or a cluster. Notice that the execution time is quiet high so this second option is highly recomended. This are the default settings for both files:

```bash
seq='WT' #protein name (WT or Shuffle)
temp=320 #Temperature
cutoff=4.0 #Cutoff distance
wind=6 #number of windows desired for WESTPA
n=100 #number of chains 
L=300 #cubic box size
steps=4000000 #number of steps of the simulation to reach equilibrium

```

Since the code must generate the initial configurations and then take them to equilibrium, this may take a while. The code also generates a `.dcd` file so you can adapt the number of steps to reach equilibrium by analysing manually the trajectory. The code generates as many initial configurations as windows (`wind`) demanded, each configuration starting with a certain number of chains already in a cluster (close to each other) but the number of chains that form the desired cluster may differ from the expected.

### Local Machine

To execute the program in your machine use the `run_ini.sh` file. This file will automatically generate the necessary directories for the code to work if they do not exist. To execute the program you should use:


```Shell
./run_ini.sh
```

### Cluster

If you want to execute the code in a cluster (i.e. powered by Slurm) then use something similar to `job.srun`. 


```Shell
sbatch job.srun
```


## RESULTS

The program will return a directory named `/config` containing the initial configurations before and after equilibrium (i.e. `top_0.50.pdb`, `top_eq\_0.50.pdb` for a system with 50\% of the chains forming the cluster). As previously mentioned, the code also generates a `.dcd` with a 100 frames from the initial position to equilibrium (following the previous example, `traj_0.50.dcd`). Use the trajectory file to set the number of steps needed by your system to reach equilibrium.

*NOTE:* _in order to visualize if the chains have reached equilibrium we used [VMD] viewer using the beta mode from graphics representation. Since we wanted a different color for each aminoacid from top to bottom, the b-factor of every particle of the chain was changed in order to do it. This directory also includes the file `b_fact.py` which changes the b factors of any `.pdb file that you may have generated from other codes._

[VMD]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

# BSTATES REPOSITORY

Following with the idea of generating different configurations in order to run the _WESTPA_ simulation, since westpa requires a `/bstates` directory with all the initial states, this can be generated using the program `bstates.py`.

## REQUIREMENTS

If you have executed on the same machine the codes to generate the `/config` directory, you will be able to run the code in order to generate the `/bstates` directory. If not, the requirements are the same as in the [previous case] plus executing the `analyse_ini.py` program in order to initialize its functions:


```Shell
python analyse_ini.py
```

Notice that it is needed that the `/config` directory exists. If not, please generate it before running `bstates.py`.

## EXECUTION

Once the `/config` directory has been created, you can use the following command to generate the `\bstates` directory setting the proper values for the variables `windows` (integer) and `seq` (string) (default are set to 6 and 'WT' respectively). After that you can run the script in the `/INITIAL_CONFIG` directory using a shell:

```Shell
python bstates.py --seq 'WT' --windows 6
```

## RESULTS

The generated directory contains all necessary files for _WESTPA_ to read and none of them should be removed! The `pcoord` used is the size of the cluster (radius) of each configuration. 

NOTE: _If you run again `bstates.py`, the previous `/bstates`directory will be permanently removed. The `bstates.txt` file is generated taking into account that each configuration has the same probability but this can be changed manually by modifying the second column of each initial configuration._

[previous case]: https://github.com/Albert2424/TFM/blob/main/INITIAL_CONFIG/README.md#requirements
