# TFM
Contains all the codes used to perform the necessary simulations for the Master Final Work.

## GENERAL GUIDELINES

### ADD A PROTEIN

For all codes used in this work, if a protein different than _WT_ and _shuffle_ must be added, this are the steps to follow:

1. Open `analyse.py` and add your protein to the proteins DataFrame (if you want to generate new initial configurations for WESTPA, add your protein to `analyse\_ini.py` too).
2. A new pickle from the new group of proteins will be generated if any of the `analyse` programs are used (or imported). If not, generate the new pickle in order to include the new protein. This can be done in example by running lines 48-51 on the `analyse\_ini.py` and executing the cell. 
3. Add the corresponding `residues.csv` of your protein or modify the existing file by adding your data.
