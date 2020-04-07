#!/bin/bash

source ~/.bashrc
goconda && conda activate py3.7.4

set -euxo pipefail

# python langevin.py  traj.xtc 1.0
# echo "0" | gmx rotmat -f aniso.xtc -s aniso.gro
# python rotdiff.py rotmat.xvg msd_angle.dat
  
i="" 
for j in {1..2}
do
    python langevin.py sig-"$i"-rep-"$j".xtc
    echo "0" | gmx rotmat -f sig-"$i"-rep-"$j".xtc -s aniso.gro
    python rotdiff.py rotmat.xvg msd-"$i"-rep-"$j".dat
done


