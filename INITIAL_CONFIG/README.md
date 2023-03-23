## INITIAL CONFIGURATION

In order to start the dynamics from different starting conditions, some specific initial distributions are generated. The aim of generating distributions with different aggregate sizes is to optimize the WESTPA simulation of the dynamics of aggregation. 

To generate the initial configurations we need three codes: _init\_config.py_, which provides a first spatial distribution with stiff peptide distribution and _simulate\_ini.py_ and _analyse\_ini.py_ which are meant to bring the configuration to equilibrium while maintaining the chains that are far away from the cluster fixed. Both _simulate\_ini.py_ and _analyse\_ini.py_ are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper _Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range_. See the code in **[here]**.



[here]: https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code


### REQUIREMENTS

The code has been created using _Python 3.10.9_. Other important packages needed are:

* OpenMM ([8.0.0])
* MDAnalysis ([2.4.2]) 
* MDtraj ([1.9.7])

[8.0.0]: http://docs.openmm.org/7.0.0/userguide/application.html
[2.4.2]: https://www.mdanalysis.org/pages/installation_quick_start/
[1.9.7]: https://www.mdtraj.org/1.9.7/installation.html

All packages may be installed using [Miniconda] using the `pip` command of the conda prompt.
`
[Miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

###EXECUTION

To execute the program in your machine use the _run\_ini.sh_ file. If you want to execute the code in a cluster (i.e. powered by Slurm) then use something similar to _job.srun_. This are the default settings for the both files:

```bash
name='WT' #protein name (WT or Shuffle)
temp=300 #Temperature
cutoff=4.0 #Cutoff distance
wind=6 #number of windows desired for WESTPA
n=100 #number of chains 
L=300 #cubic box size
steps=4000000 #number of steps of the simulation to reach equilibrium

```



