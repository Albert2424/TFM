### Initial Configurations

In order to start the dynamics from different starting conditions, some specific initial distributions are generated. The aim of generating distributions with different aggregate sizes is to optimize the WESTPA simulation of the dynamics of aggregation. 

To generate the initial configurations we need three codes: _init\_config.py_, which provides a first spatial distribution with stiff peptide distribution and _simulate\_ini.py_ and _analyse\_ini.py_ which are meant to bring the configuration to equilibrium while maintaining the chains that are far away from the cluster fixed. Both _simulate\_ini.py_ and _analyse\_ini.py_ are adaptations of codes made by Giulio Tesei and Kresten Lindorff-Larsen on their paper [Improved Predictions of Phase Behaviour of IDPs by Tuning the Interaction Range](https://github.com/KULL-Centre/papers/tree/main/2022/CG-cutoffs-Tesei-et-al/MC/code.
