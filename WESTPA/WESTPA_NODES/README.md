# WESTPA SIMULATION OF PROTEIN AGREGATION (USING NODES)

This repository contains all the necessary files to perform the _[WESTPA] simulation of protein agregation_ (WSPA).

[WESTPA]: https://github.com/westpa/westpa

## TABLE OF CONTENTS
1. [ Requirements ](#1-req)
2. [ Contains ](#2-contains)

	2.1 [ `/bstats` ](#2.1-bs)
	
	2.2 [ `/common_files` ](#2.2-cf)
	
	2.3 [ `/westpa_scripts` ](#2.3-ws)
	
	2.4 [ `env.sh` and `ini.sh` ](#2.4-eir)
	
	2.5 [ `west.cfg` ](#2.5-w)
	
	2.6 [ `runwe.slurm` ](#2.6-r)
	
	2.7 [ `input.dat` and `set_input.sh` ](#2.7-is)
	
3. [ Execution ](#3-e)

	3.1 [ Choosing Parameters ](#3.1-cp)

	3.2 [ Cluster ](#3.1-c)
	
	3.3 [ Local Machine ](#3.2-lm)

4. [ Results ](#4-r)
	
	4.1 [ Plots using wedap ](#4.1-wedap)
	
	4.2 [ Plots using w_pdist and plothist ](#4.2-w_pdist)


<a name="1-req"></a>
## Requirements

Since _WESTPA_ will be used, follow the specific [requirements] of the library. In addition, for the simulation codes and the posterior analysis other important packages needed are:

* OpenMM ([8.0.0])
* MDAnalysis ([2.4.2]) 
* MDtraj ([1.9.7])
* wedap ([0.0.6])

[8.0.0]: http://docs.openmm.org/7.0.0/userguide/application.html
[2.4.2]: https://www.mdanalysis.org/pages/installation_quick_start/
[1.9.7]: https://www.mdtraj.org/1.9.7/installation.html
[0.0.6]: https://pypi.org/project/wedap/

All packages may be installed using [Miniconda] with the `pip` command of the conda prompt.

[Miniconda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Make sure that all this packages are included in your _WESTPA_ enviorenment. Strictly speaking it is not necessary to have multiple nodes to run the simulation using this code, however if you only have one node it is better to use the code in [WESTPA\_GPU]. All nodes must have the same number of GPUs or, at least, the maximum number of GPUs used must be maximum as high as the number of GPUs of the node with less GPUs.

To change the number of nodes used modify the `runwe.slurm` file:

```Shell
#SBATCH -N 2
```

[requirements]: https://github.com/westpa/westpa#requirements
[WESTPA\_GPU]: https://github.com/Albert2424/TFM/edit/main/WESTPA/WESTPA_GPU

<a name="2-contains"></a>
## Contains
<a name="2.1-bs"></a>
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


<br></br>
<a name="2.2-cf"></a>
- **`/common_files`:** This directory contains the necessary _OpenMM_ files required to perform the simmulation and also the `bstate.pdb` file which contains the topology of the system. Eventhough it seems like it is not being used in the _openMM_ simulation, the `bstate.pdb` is required for _WESTPA_ to know the topology of the system. Basically you can use any `.pdb` file which contains a frame of your system. The default file is equivalent to the equilibrated configuration with no cluster formation. To modify this directory follow these steps:

1. Add all necessary files to perform the simulation. In our case they are
`simulate.py`, `analyse.py` and `residues.csv`. Notice that the `protein.pkl` will be generated automatically in your machine to avoid errors due to binary files.
2. Change the `bstate.pdb` to a _pdb_ that contains the topology of your system.
3. Notice that the tau value is defined at the simulation files and it corresponds to the time that each simulation is ran.
4. Add the file that will calculate the pcoord of every simulation. In this case this is done by `cluster.py` wich generates a file named `dist.dat` containing the pcoord.

Both `simulate.py` and `analyse.py` are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper _Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range_. See the original code in **[here]**.

[here]: https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code

<br></br>
<a name="2.3-ws"></a>
- **`/westpa_scripts`:** This directory contains all necessary files to perform the _WESTPA_ simulation. Most of them are default files that must not be modified. However, to make _WESTPA_ perform the desired simulation we need to modify the `runseg.sh` file that controlls how each segment is simulated. In order to do it you must follow these steps:

1. Change every _simulate.py_ for your simulation code. Notice that when the simulation code is called some arguments are added. If your code does not need this arguments you must erase them. On the top of the `runseg.sh` file the necessary variables needed in order to perform the simulation are defined.   
2. Change _cluster.py_ for your code that gets the pcoord. As in the previous case some arguments where added.
3. You need to deliver to _WESTPA_ all the necessary files to perform the simulation. you can do this adding (or eraseing) commands like the following: `sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/file_to_be_added > file_to_be_added`
4. If the name of your basis states (and generated checkpoints) is not _top.pdb_ then change all _top.pdb_ for the desired name. 

<br></br>
<a name="2.4-eir"></a>
- **`env.sh`** and **`ini.sh`:** This are the needed files to ensure that the _WESTPA_enviorenment is set correctly and to initialize and run the simulation. Since in this case we do not have a target state, the lines containing _TSTATE\_ARGS_ on the _w\_init_ command have been commented. If you need a target state, then uncomment these lines and specify your target state on the `tstate.file` (in the default simulation this file is useless).  

<br></br>
<a name="2.5-w"></a>
- **`west.cfg`:** This is the file containing the configuration parameters for the westpa simulation. For more information look [here]. In this case [addaptative binning] (MAB) was used. The MAB needs a starting bin set. You can use the [bstates.py] program to obtain a bin set suggestion or you can modify it manually. The number of starting bins is not fixed so you can refine the bin set for your specific simulation. In the default file, the boundaries are `[0.0,inf]`. 

You should also consider changing the dimension and number of data points of the pcoord (_pcoord\_dim_ and _pcoord\_len_ respectively). 



[bstates.py]: https://github.com/Albert2424/TFM/blob/main/INITIAL_CONFIG/bstates.py
[here]: https://westpa.github.io/westpa/users_guide/west/setup.html#configuration-file
[addaptative binning]: https://westpa.readthedocs.io/en/latest/users_guide/west/setup.html#recursivebinmapper

<br></br>
<a name="2.6-r"></a>
- **`runwe.slurm`:** Example file to run your simulation on  a cluster powered by slurm. Notice that the name of the job will always be the name of the repository containing the simulation.

<br></br>
<a name="2.7-is"></a>
- **`input.dat`** and **`set_input.sh`:** Files containig the most relevant input changes that may be modified for each simulation. 
<a name="3-e"></a>
## Execution
<a name="3.1-cp"></a>
### Choosing parameters

Before starting the simulation you may need to adjust some parameters such as the radius at wich a particle is included in a cluster or not. This can be done using the `rc_graph.py` program in the `/extra` directory with a bunch of configurations. The default program is made to test the radius with the initial configurations generated with the [ `bstates.py` ] program which are stored in the `/config` directory. Execute the program in the `/extra` directory to obtain a graph using:

```Shell
python rc_graph.py --seq 'WT' --windows 15 --L 300. --n_chains 100
```
Where `seq` is the name of your protein, `windows` is the number of configurations you have, `L` is the length of the simulation box and `n_chains` is the number of chains of the simulation. 

[ `bstates.py` ]: https://github.com/Albert2424/TFM/blob/main/INITIAL_CONFIG/bstates.py

<a name="3.1-c"></a>
### Cluster

Prior to starting the siulation you may want to change the parameters of the simulation. This can be done by modifying the variables of the `input.dat` file. The input file is divided in 3 section: one for the simulation parameters, one for the westpa configuration and one for the cluster configuration. This are the default settings:


```diff
############### SIMULATION INPUT ############################

seq=⏩'WT'⏪       #name of the studied protein (str)
temp=⏩320⏪       #temperature (int)
cutoff=⏩4.0⏪     #cutoff (float)
tau=⏩4000⏪       #steps of the simulation (tau in westpa) (int)
n_chains=⏩100⏪   #number of chains of the system (int)
L=⏩300.⏪         #length of the box (float)
rc=⏩30.⏪         #minimum radius to consider a chain is inside a cluster (float)

############### WESTPA CONFIGURATION ########################

#number of required westpa iterations
    max_total_iterations: ⏩100⏪
#same (or higher) as slurm --time
    max_run_wallclock:    ⏩72:00:00⏪
#name of the westpa .h5 file
    west_data_file: ⏩west.h5⏪
#number of adaptative bins used:
            nbins: [⏩30⏪]

############### SLURM #######################################

#SBATCH --time=⏩08:00:00⏪
```

Notice that only the values between arrows (`⏩⏪`) may be changed for the code to work properly. If for any reason you need to modify the `ẁest.cfg` file (or any of the other two), check if the line number of the file has changed. If it is the case, modify the `set_input.sh` file in order for the lines to fit as suggested in this example:

```
sed -i '⏩31⏪s/.*/'"$(sed -n '12p' input.dat)"'/' west.cfg   <-- set_input.sh

31.     max_total_iterations: 100                            <-- west.cfg 

``` 

To execute the program in a cluster (i.e. powered by Slurm) use something similar to `runwe.slurm` and the `Makefile` using:


```
make run
``` 

or 

```
make
```

If you want to cancel you job you may use `make cancel`. In addition you can clean the results of a previous simulation using `make clean`. To look up for errors in the process look at the different `.log` files generated, the `.err` file or in the `/seg_log` directory that is generated with the _logs_ of every segment. If you are debugging your program consider seting to 1 the debug options in the `west.cfg` file. 

<a name="3.2-lm"></a>
### Local Machine

Eventhough it is not as efficient as the use of a cluster, if you want to do some testing using your machine you may want to use:

```shell
./init.sh
```
Followed by:

```shell
./run.sh
```

If you need to update the input, use (before the previous steps):

```shell
./set_input.sh
```

Notice that you still need the `runwe.slurm` file or you can erase the last line of the `set_input.sh` file.
<a name="4-r"></a>
## Results

Once the simulation finishes, a file named `west.h5` (or the name you have set in the `input.dat` file) is generated. This file has all the information of the simulation and can be analysed using the _westpa_ package itself. However, another option is to use the [ `wedap` package ] which basically makes it easier to make the plots.

_NOTE: the following sections consider that the name of the _westpa_ file is the default `west.h5`. If it is not your case, you may change this name in the `Makefile`._

[ `wedap` package ]: https://pypi.org/project/wedap/

<a name="4.1-wedap"></a>
### Plots using wedap 

The fastest way to plot both the average plot and evolution plot is to use the Makefile option in the same directory containing the `west.h5` file:

```Shell
make analyse
``` 
This will generate `hist_evo.pdf` and `hist_av.pdf` which are the 1D energy evolution for every iteration including the value of the pcoord and the average energy for every pcoord respectively. If you need different options use:

```Shell
wedap --help
```
To display all the available options.  

<a name="4.2-w_pdist"></a>
### Plots using w_pdist and plothist 

Another option is to use the tools provided by _westpa_ itself. You can rapidly plot what is needed if you use:


```Shell
w_pdist
```

And then use the generated `pdist.h5` file with plothist just like:


```Shell
plothist evolution pdist.h5 -o hist_evo.pdf
```

```Shell
plothist average pdist.h5 -o hist_av.pdf
```

Use the `--help` option in both cases to obtain more information or look at the [ westpa manual ].

[ westpa manual ]: https://westpa.readthedocs.io/_/downloads/en/stable/pdf/
