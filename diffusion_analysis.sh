#!/bin/bash

source ~/.bashrc
goconda && conda activate py3.7.4

SIG_TRANS=2
SIG_ROT=0.1

set -euxo pipefail

# run N replicas for better statistics
for i in {1..1}
do
    # perform the simulations
    python langevin.py traj_comp.$i.xtc $SIG_TRANS $SIG_ROT

    # rotation matrix of the original trajectory
    echo "1" | gmx rotmat -f traj_comp.$i.xtc \
                          -s aniso.gro \
                          -o rotmat.$i.xvg

#   # rotational MSD
    python rotdiff.py rotmat.$i.xvg msr.$i.dat

#   # Anisotropy
    python diffusion_tensor.py aniso.gro traj_comp.$i.xtc rotmat.$i.xvg Dtensor.$i.dat
    python orient.py Dtensor.$i.dat Dtensor-DIAG.$i.dat

done
