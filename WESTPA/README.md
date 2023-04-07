# WESTPA SIMULATION OF PEPTIDE AGREGATION

This repository contains all the necessary files to perform the _[WESTPA] simulation of peptide agregation_ (WSPA).

[WESTPA]: https://github.com/westpa/westpa
## CONTAINS

`/bstates`: Repository containing the `bstates.txt` and `pcoord.init` files and the repositories containing every initial configuration (in `.pdb` format) and a `pcoord.init` file with the reaction coordinate of said configuration. Here's a brief explanation on how to create manually the `/bstates` directory: 
      1. Add all the basis states you may need (initial configurations) in                different directories.
      2. Add a `pcoord.init` file to each directory containning the reaction              coordinate (pcoord) of the corresponding initial configuration.
      3. Open `bstates.txt` and add all the directories that are part of the              simulation:
            ID      probability     filename
            0       0.17            dir1
            1       0.17            dir2
            .       .               .
            .       .               .
            .       .               .
         
         The probability is automatically renormalized to 1 so it does not                matter if it does not add up to 1.
         
       4. Open `pcoord.init` file from the `/bstates` directory and add the               pcoord of every bstate.
