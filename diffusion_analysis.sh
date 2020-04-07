#!/bin/bash

source ~/.bashrc
goconda && conda activate py3.7.4

set -euxo pipefail

# python langevin.py  traj.xtc 1.0
# echo "0" | gmx rotmat -f aniso.xtc -s aniso.gro
# python rotdiff.py rotmat.xvg msd_angle.dat
  
# Translation
for i in {1..5}
do
    # Rotation
    for j in {1..5}
    do
        python langevin.py trans-$i-rot-0.$j.xtc $i 0.$j
        echo "0" | gmx msd -xvg none -f trans-$i-rot-0.$j.xtc -s aniso.gro -lateral z -o msdt_trans-$i-rot-0.$j.xvg
        echo "0" | gmx rotmat -f trans-$i-rot-0.$j.xtc -s aniso.gro
        python rotdiff.py rotmat.xvg msdr-$i-rot-0.$j.dat
    done
done
