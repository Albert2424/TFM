## INITIAL CONFIGURATION

In order to start the dynamics from different starting conditions, some specific initial distributions are generated. The aim of generating distributions with different aggregate sizes is to optimize the WESTPA simulation of the dynamics of aggregation. 

To generate the initial configurations we need three codes: `init_config.py`, which provides a first spatial distribution with stiff peptide distribution and `simulate_ini.py` and `analyse_ini.py` which are meant to bring the configuration to equilibrium while maintaining the chains that are far away from the cluster fixed. Both `simulate_ini.py` and `analyse_ini.py` are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper _Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range_. See the original code in **[here]**.



[here]: https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code


### REQUIREMENTS

The code has been created using _Python 3.10.9_. Other important packages needed are:

* OpenMM ([8.0.0])
* MDAnalysis ([2.4.2]) 
* MDtraj ([1.9.7])

[8.0.0]: http://docs.openmm.org/7.0.0/userguide/application.html
[2.4.2]: https://www.mdanalysis.org/pages/installation_quick_start/
[1.9.7]: https://www.mdtraj.org/1.9.7/installation.html

All packages may be installed using [Miniconda] with the `pip` command of the conda prompt.

[Miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html


### EXECUTION

To execute the program in your machine use the `run_ini.sh` file. If you want to execute the code in a cluster (i.e. powered by Slurm) then use something similar to `job.srun`. This are the default settings for the both files:

```bash
seq='WT' #protein name (WT or Shuffle)
temp=300 #Temperature
cutoff=4.0 #Cutoff distance
wind=6 #number of windows desired for WESTPA
n=100 #number of chains 
L=300 #cubic box size
steps=4000000 #number of steps of the simulation to reach equilibrium

```

Since the code must generate the initial configurations and then take them to equilibrium, this may take a while. The code also generates a `.dcd` file so you can adapt the number of steps to reach equilibrium. The code generates as many initial configurations as windows (`wind`) demanded, each configuration with a certain number of chains already in a cluster.


### RESULTS

The program will return a directory named `/config` containing the initial configurations before and after equilibrium (i.e. `top_0.50.pdb`, `top_eq\_0.50.pdb` for a system with 50\% of the chains forming the cluster). As previously mentioned, the code also generates a `.dcd` with a 100 frames from the initial position to equilibrium (following the previous example, `traj_0.50.dcd`). Use the trajectory file to set the number of steps needed by your system to reach equilibrium.

*NOTE:* _in order to visualize if the chains have reached equilibrium we used [VMD] viewer using the beta mode from graphics representation. Since we wanted a different color for each aminoacid from top to bottom, the b-factor of every particle of the chain was changed in order to do it. This directory also includes the file `b_fact.py` which changes the b factors of any `.pdb file that you may have generated from other codes._

[VMD]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD





