# WESTPA SIMULATION OF PEPTIDE AGREGATION

This repository contains all the necessary files to perform the _[WESTPA] simulation of peptide agregation_ (WSPA).

[WESTPA]: https://github.com/westpa/westpa

## REQUIREMENTS

Since _WESTPA_ will be used, follow the specific [requirements] of the library. In addition, for the simulation codes other important packages needed are:

* OpenMM ([8.0.0])
* MDAnalysis ([2.4.2]) 
* MDtraj ([1.9.7])

[8.0.0]: http://docs.openmm.org/7.0.0/userguide/application.html
[2.4.2]: https://www.mdanalysis.org/pages/installation_quick_start/
[1.9.7]: https://www.mdtraj.org/1.9.7/installation.html

All packages may be installed using [Miniconda] with the `pip` command of the conda prompt.

[Miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Make sure that all this packages are included in your _WESTPA_ enviorenment.

[requirements]: https://github.com/westpa/westpa#requirements

## CONTAINS

- **`/bstates`:** Repository containing the `bstates.txt` and `pcoord.init` files and the repositories containing every initial configuration (in `.pdb` format) and a `pcoord.init` file with the reaction coordinate of said configuration. Here's a brief explanation on how to create manually the `/bstates` directory: 
 1. Add all the basis states you may need (initial configurations) in different directories. All the basis states must have the same name, in this case `top.pdb` (if you want to use a different name, you must change other files for your simulation to work).
 2. Add a `pcoord.init` file to each directory containning the reaction coordinate (pcoord) of the corresponding initial configuration.
 3. Open `bstates.txt` and add all the directories that are part of the simulation:
  ```
   ID      probability     filename
   0       0.17            dir1
   1       0.17            dir2
   .       .               .
   .       .               .
   .       .               .
  ```
  The probability is automatically renormalized to 1 so it does not matter if it does not add up to 1.
         
4. Open `pcoord.init` file from the `/bstates` directory and add the pcoord of every basis state.

Notice that the `/bstates` directory can be created using the program [bstates.py] on the `/INITIAL_CONFIG` directory.

[bstates.py]: https://github.com/Albert2424/TFM/blob/main/INITIAL_CONFIG/bstates.py


- **`/common_files`:** This directory contains the necessary _OpenMM_ files required to perform the simmulation and also the `bstate.pdb` file which contains the topology of the system. Eventhough it seems like it is not being used in the _openMM_ simulation, the `bstate.pdb` is required for _WESTPA_ to know the topology of the system. Basically you can use any `.pdb` file which contains a frame of your system. The default file is equivalent to the equilibrated configuration with no cluster formation. To modify this directory follow these steps:

1. Add all necessary files to perform the simulation. In our case they are
`simulate.py`, `analyse.py` and `residues.csv`. Notice that the `protein.pkl` will be generated automatically in your machine to avoid errors due to binary files.
2. Change the `bstate.pdb` to a pdb that contains the topology of your system.
3. Notice that the tau value is defined at the simulation files and it corresponds to the time that each simulation is ran.
4. Add the file that will calculate the pcoord of every simulation. In this case this is done by `cluster.py` wich generates a file named `dist.dat` containing the pcoord.

Both `simulate.py` and `analyse.py` are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper _Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range_. See the original code in **[here]**.

[here]: https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code

- **`/westpa_scripts`:** This directory contains all necessary files to perform the _WESTPA_ simulation. Most of them are default files that must not be modified. However, to make _WESTPA_ perform the desired simulation we need to modify the `runseg.sh` file that controlls how each segment is simulated. In order to do it you must follow these steps:

1. Change every _simulate.py_ for your simulation code. Notice that when the simulation code is called some arguments are added. If your code does not need this arguments you must erase them. On the top of the `runseg.sh` file the necessary variables needed in order to perform the simulation are defined.   
2. Change _cluster.py_ for your code that gets the pcoord. As in the previous case some arguments where added.
3. You need to deliver to _WESTPA_ all the necessary files to perform the simulation. you can do this adding (or eraseing) commands like the following: `sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/file_to_be_added > file_to_be_added`
4. If the name of your basis states (and generated checkpoints) is not _top.pdb_ then change all _top.pdb_ for the desired name. 


- **`env.sh`**,**`ini.sh`** and **`run.sh`:** This are the needed files to ensure that the _WESTPA_enviorenment is set correctly and to initialize and run the simulation. Since in this case we do not have a target state, the lines containing _TSTATE\_ARGS_ on the _w\_init_ command have been commented. If you need a target state, then uncomment these lines and specify your target state on the `tstate.file` (in the default simulation this file is useless).  

- **`west.cfg`:** This is the file containing the configuration parameters for the westpa simulation. For more information look [here]. In this case [addaptative binning] (MAB) was used. The MAB needs a starting bin set. You can use the [bstates.py] program to obtain a bin set suggestion or you can modify it manually. The number of starting bins is not fixed so you can refine the bin set for your specific simulation. In the default file, the boundaries are `[0.0,0.01,5.24,10.47,15.71,20.95,26.18,inf]`. 

You should also consider changing the dimension and number of data points of the pcoord (_pcoord\_dim_ and _pcoord\_len_ respectively). 



[bstates.py]: https://github.com/Albert2424/TFM/blob/main/INITIAL_CONFIG/bstates.py
[here]: https://westpa.github.io/westpa/users_guide/west/setup.html#configuration-file
[addaptative binning]: https://westpa.readthedocs.io/en/latest/users_guide/west/setup.html#recursivebinmapper

- **`runwe.slurm`:** Example file to run your simulation on  a cluster powered by slurm.

## EXECUTION

To execute the program in a cluster (i.e. powered by Slurm) then use something similar to `runwe.slurm` and the `Makefile` using `make run`. 

Prior to starting the siulation you may want to change the parameters of the simulation. This can be done by modifying the variables of the `runseg.sh` file in the `/westpa_scripts` directory. This are the default settings:

```bash
name='WT'      #name of the studied protein (str)
temp=320       #temperature (int)
cutoff=4.0     #cutoff (float)
tau=400        #steps of the simulation (tau in westpa) (int)
n_chains=100   #number of chains of the system (int)
```
To look up for errors in the process look at the different `.log` files generated, the `.err` file or in the `/seg_log` directory that is generated with the _logs_ of every segment.

## RESULTS



