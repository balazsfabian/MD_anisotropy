import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
import sys

# np.random.seed(seed=0)
##################################
# SIMULATING THE LANGEVIN EQUATION
#    FOR ANISOTROPIC PARTICLES
#
# REF: Ribault, Triller, Sekimoto
#  Phys. Rev. E 75, 021112 (2007)
##################################


def Mt (theta, eta):
    """ Matrix determining the orientation of displacements.
        When applied to [Dx, Dy].T, it produces displacements
        that have the same norm as [Dx, Dy].T
        theta : orientation in radians
        eta   : anisotropy
    """
    ut = np.array([np.cos(theta), np.sin(theta)])
    mat = np.outer(ut,ut)

    M = np.sqrt(1+eta) * mat + np.sqrt(1-eta)*(np.eye(2)-mat) 

    return M


def update_step(x,a, eta, x_var=1.0, a_var=0.1):
    """ Update anisotropic Langevin position and angle
        x_var : variance of the step distribution
        y_var : variance of the angle distribution
    """

    # Generating the random variables
    # with given variances
    W = np.random.normal(scale=np.sqrt(a_var))
    B = np.random.normal(scale=np.sqrt(x_var),size=2)

    # Updating postions and angle
    # Note: this way, Mt is not stored and
    # needs to be recalculated at every step.
    M_old = Mt(a, eta)
    a_new = a + W
    M_new = Mt(a_new, eta)
    x_new = x + 0.5 * np.dot(M_new + M_old,B)

    return x_new, a_new


def create_diatom(x, l, theta):
    """ Create 3D positions from 2D center of mass, length and angle
    """
    c = np.cos(theta)
    s = np.sin(theta)
    dx  = np.array([c, s]) * l

    pos = np.zeros((3,3))
    pos[0,0:2] = x + dx
    pos[1,0:2] = x - dx
    pos[2,0:2] = x
    pos[2,2]   = 1.0     # elevation of the COM
    return pos

# Parameters of the trajectory
#   On diffusion coefficients:
#  * MSD_rot   = a_var * t = 2D * t
#  * MSD_trans = 2 * x_var * t = 4D * t
#   Moreover, gmx output is in nm^2
#   instead of angstrom^2
#  WARNING: setting a_var too large
#   results in undersampling of the
#   rotational motion and it causes
#   the underestimation of MSD_rot.
#   (e.g.: 2pi => 0)
#
eta     =       0  # degree of anisotropy
x_var   =       1  # variance of the step distribution (in angstrom^2)
a_var   =     0.1  # variance of the angle distribution (in radian^2),
                   # a_var = \sigma_{angle}^2 = 2D
n_steps =  100000  # steps
angle   =       0  # initial angle in radians
com     = [0.,0.]  # initial position
length  =       1  # length of the particle in angstroms
l_box   =     100  # box dimensions in angstrom

# Create MDA Universe from scratch:
atoms    = 3
residues = 1
u = mda.Universe.empty(atoms,
                       residues,
                       atom_resindex=[0, 0, 0],
                       trajectory=True)

u.add_TopologyAttr('names')
u.add_TopologyAttr('resids')
u.add_TopologyAttr('resnames')
u.atoms.names = ['C1','C2','N']
u.residues.resids = "1"
u.residues.resnames = "ANI"
coordinates = np.empty((n_steps, u.atoms.n_atoms,3))


# Generate trajectory
for i in tqdm(range(len(coordinates))):

    # create the positions of a diatomic molecule
    # with given length and store them
    pos = create_diatom(com, length, angle)
    coordinates[i] = pos

    # update the positions
    com, angle = update_step(com, angle, eta, x_var=float(sys.argv[2]), a_var=float(sys.argv[3]))

# load the coordinates into the universe
# and write them in a GRO and XTC format
u.load_new(coordinates, order='fac')
u.dimensions = np.array([l_box,l_box, l_box, 90, 90, 90])
u.atoms.write('aniso.gro')
with mda.Writer(sys.argv[1], u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        u.dimensions = np.array([l_box,l_box, l_box, 90, 90, 90])
        W.write(u.atoms)

