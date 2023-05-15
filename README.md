# TFM
Contains all the codes used to perform the necessary simulations for the Master Final Work.

## General Guidelines

### Add a Protein

For all codes used in this work, if a protein different than _WT_ and _shuffle_ must be added, this are the steps to follow:

1. Open `analyse.py` and add your protein to the proteins DataFrame (if you want to generate new initial configurations for WESTPA, add your protein to `analyse_ini.py` too).
2. A new pickle from the new group of proteins will be generated if any of the `analyse` programs are used (or imported). If not, generate the new pickle in order to include the new protein. This can be done in example by running lines 48-51 on the `analyse_ini.py` and executing the cell. 
3. If you want to use different parameters for the residues, add the corresponding `residues.csv` to the directory or modify the existing file by adding your data.


### Contents

1. `/INITIAL_CONFIG`: Directory containning all the necessary files to generate the initial configurations of the system to use in the WESTPA enviorement.
2. `/WESTPA`: Directory containning all the necessary files and directories to perform a _WESTPA_ simulation of peptide agregation.  

