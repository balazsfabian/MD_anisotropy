import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import MDAnalysis as mda
from tqdm import tqdm


def TDT(pos, ori, max_dt=400):

    # pos and ori must have the same length
    l_traj = len(pos)
    assert (l_traj == len(ori)), "ERROR: mismatch in trajectory and rotmat length"

    # project pos and ori onto the XY plane
    pos = pos[:,:2]
    ori = ori[:,:2,:2]

    # max_dt cannot be larger than
    # the length of the trajectory
    max_dt = min(l_traj, max_dt+1)

    # array for the diffusion tensor
    # and the normalization
    MSD_ij = np.zeros((max_dt,2,2))
    norm   = np.zeros((max_dt,2,2))

    # initial time t0 of the average
    for t0 in tqdm(range(l_traj-1)):

        # determine the maximum length
        # of trajectory from t0 onward
        max_dt = min(max_dt, len(pos)-t0)

        # rotation needed to align the
        # first frame of the trajectory
        # with the reference.
        rot  = ori[t0]
 
        for dt in range(1, max_dt):

            # rotate the displacement into
            # the reference orientation
            dx = np.dot(rot, pos[dt+t0]-pos[t0])
 
            # calculate MSD and increment
            # the normalization
            MSD_ij[dt] += np.outer(dx,dx)
            norm[dt]   += 1

    # increment the 0 bin of the normalization,
    # as the corresponding MSD_ij[0] is also zero,
    # and apply the normalization
    norm[0] += 1
    MSD_ij   = np.divide(MSD_ij,norm)
 
    return  MSD_ij


# load the universe, trajectory and
# create a new auxiliary attribute
# containing the rotational data
u = mda.Universe(sys.argv[1], sys.argv[2])
u.trajectory.add_auxiliary('orientation', sys.argv[3])

# array for positions and orientations
frames = u.trajectory.n_frames
pos = np.zeros((frames,3))
ori = np.zeros((frames,3,3))
dt  = u.trajectory.dt
sel = u.select_atoms('resname Ar4 or resname ANI')
print (sel)

# collect the relevant data from the trajectory
for i,ts in enumerate(u.trajectory):

    # load the center of mass of the particle 
    # into an array. For this to work, the
    # trajectory must only contain the particle
    # or its center of mass
    pos[i] = sel.center_of_geometry()
#   pos[i] = sel[0].position

    # load the rotation matrix as a row vector
    # and reshape it into a 3D rotation matrix.
    # This matrix represents the rotation needed
    # to align the configuration with the reference.
    raw_rot = ts.aux.orientation[1:]
    ori[i]  = raw_rot.reshape((3,3))


# calculate the time-dependent diffusion tensor
# at a fixed given orientation. This orientation
# is the one specified by the rotmat.xvg
# Then, reshae the array and save it
MSD_ij = TDT(pos, ori, max_dt=10000)
MSD_ij = MSD_ij.reshape((len(MSD_ij),4))

# set up the lagtime array, reshape it and
# attach it as a first column to the data
l_data  = len(MSD_ij)
lagtime = np.linspace(0, l_data-1, l_data) * dt
lagtime = lagtime.reshape((len(lagtime),1))
result  = np.hstack((lagtime,MSD_ij))

# write the final array
np.savetxt(sys.argv[4], result, header = "time(ps) MSD_xx(angstrom^2) MSD_xy(angstrom^2) MSD_yx(angstrom^2) MDS_yy(angstrom^2)")
