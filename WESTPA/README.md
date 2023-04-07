# WESTPA SIMULATION OF PEPTIDE AGREGATION

This repository contains all the necessary files to perform the _[WESTPA] simulation of peptide agregation_ (WSPA).

[WESTPA]: https://github.com/westpa/westpa
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

- **`/westpa_scripts`:**

