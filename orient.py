import sys
import numpy as np
from scipy.spatial.transform import Rotation as R

# load diffusion tensor that was
# calculated at a fixed angle.
# Assign time and diffusion tensor
# to variables.
raw_data = np.loadtxt(sys.argv[1])
time = raw_data[:,0]
MSD_ij = raw_data[:,1:]

# Reshape the row vector to a 2D matrix
MSD_ij = MSD_ij.reshape((len(MSD_ij),2,2))

# at every step (lagtime), determine the eigevalues and eigenvectors of the 2D
# diffusion tensor. Here, b[:,i] is the eigenvector corresponding to the eigenvalue a[i].
# The eigenvalues at a given lagtime are the MSD_xx and MSD_yy that maximize the
# difference between the major and minor axes of the particle. The major axis aligns with
# the laboratory frame's X, the minor with the Y. The eigenvalues are in "a". The rotation
# matrix needed to accomplish this is contained in the eigenvector matrix.
a,b = np.linalg.eig(MSD_ij)

# sort the MSD-s (that is, the eigenvalues)
# and the corresponding vectors
for i,elem in enumerate(a):
    idx = elem.argsort()[::-1]
    a[i] = a[i][idx]
    b[i] = b[i][:,idx]


# at every step (lagtime), calculate the angle from
# the eigenvector associated with the major/minor axes.
angle = np.zeros(a.shape)
for i,elem in enumerate(b):

    # w is the array of eigenvectors
    # w[:,i] is the eigenvector of a[i]
    w = b[i]
    for j in range(2):
        angle[i,j] =  np.arctan2(w[1,j], w[0,j])


# Unwrap the angles to be continuous.
# unwrap works in 180 degree complements,
# but the orientation of the eigenvectors
# is in 90 degree complements, hence the
# multiplication and division by 2.
angle[:,0] = np.unwrap(angle[:,0]*2)/2
angle[:,1] = np.unwrap(angle[:,1]*2)/2
 
np.savetxt(sys.argv[2], np.vstack((time.T, a.T,angle.T)).T, header = "time(ps) MSD_maj(angstrom^2) MSD_min(angstrom^2) angle_maj(rad) angle_min(rad)")
