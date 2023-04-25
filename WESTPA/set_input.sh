#!/bin/bash

sed -i '5s/.*/'"$(sed -n '3p' input.dat)"'/' westpa_scripts/runseg.sh
sed -i '6s/.*/'"$(sed -n '4p' input.dat)"'/' westpa_scripts/runseg.sh
sed -i '7s/.*/'"$(sed -n '5p' input.dat)"'/' westpa_scripts/runseg.sh
sed -i '8s/.*/'"$(sed -n '6p' input.dat)"'/' westpa_scripts/runseg.sh
sed -i '9s/.*/'"$(sed -n '7p' input.dat)"'/' westpa_scripts/runseg.sh

sed -i '31s/.*/'"$(sed -n '12p' input.dat)"'/' west.cfg
sed -i '32s/.*/'"$(sed -n '14p' input.dat)"'/' west.cfg
sed -i '36s/.*/'"$(sed -n '16p' input.dat)"'/' west.cfg
sed -i '23s/.*/'"$(sed -n '18p' input.dat)"'/' west.cfg

sed -i '9s/.*/'"$(sed -n '22p' input.dat)"'/' runwe.slurm
