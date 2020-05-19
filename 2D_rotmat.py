import sys
import MDAnalysis as mda
import numpy as np
from scipy.spatial.transform import Rotation as R

# GMX-STYLE ROTMAT FOR 1D RODS IN 2D SIMULATIONS

def header_lines():
    text = """# gmx rotmat-style rotmat from 2D simulation
#
# S  C  A  M  O  R  G
#
@    title "Fit matrix"
@    xaxis  label "Time (ps)"
@    yaxis  label ""
@TYPE xy
@ view 0.15, 0.15, 0.75, 0.85
@ legend on
@ legend box on
@ legend loctype view
@ legend 0.78, 0.8
@ legend length 2
@ s0 legend "xx"
@ s1 legend "xy"
@ s2 legend "xz"
@ s3 legend "yx"
@ s4 legend "yy"
@ s5 legend "yz"
@ s6 legend "zx"
@ s7 legend "zy"
@ s8 legend "zz"""
    return text

def vecs2angle(vector_1,vector_2):
    """ Angle between two vectors """

    vector_1 = vector_1[:2]
    vector_2 = vector_2[:2]
    
    angle = np.arctan2(vector_2[1], vector_2[0]) - np.arctan2(vector_1[1], vector_1[0])

    return angle


# Fixed reference
u_ref = mda.Universe(sys.argv[1])
sel = u_ref.select_atoms('name MAr1')
ref_pos = sel.positions
ref_vec = ref_pos[1] - ref_pos[0]

# Trajectory
u = mda.Universe(sys.argv[1], sys.argv[2])
sel = u.select_atoms('name MAr1')

# Array for results
n_frames = u.trajectory.n_frames
result = np.zeros((n_frames,10))

# Loop over the trajectory
for ts in u.trajectory:

    # vector to angle
    pos = sel.positions
    vec = pos[1] - pos[0]
    angle = vecs2angle(vec, ref_vec)
    
    # angle to rotmat
    rotobj = R.from_rotvec(angle*np.array([0,0,1]))
    rotmat = rotobj.as_dcm()

    # fill results
    result[ts.frame][0] = ts.time
    result[ts.frame][1:] = rotmat.flatten()

np.savetxt(sys.argv[3], result, header=header_lines(),fmt="%6d"+9*"%8.4lf",comments="")
