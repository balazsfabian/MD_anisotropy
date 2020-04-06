import numpy as np
import MDAnalysis as mda
from tqdm import tqdm
import sys

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
    #TODO: check the normalization of this matrix for various eta
    ut = np.array([np.cos(theta), np.sin(theta)])
    mat = np.outer(ut,ut)

    M = np.sqrt(1+eta) * mat + np.sqrt(1-eta)*(np.eye(2)-mat) 

    return M


def update_step(x,a, eta, x_size=1.0, x_center = 0.0, a_size=1.0, a_center= 0.0):
    """ Update anisotropic Langevin position and angle
    """

    # Generating the random variables
    W = np.random.normal()
    W = np.radians(W)
    B = np.random.normal(0,1,2)

    # Updating postions and angle
    # Note: this way, Mt is not stored and
    # needs to be recalculated at every step.
    M_old = Mt(a, eta)
    a_new = a + np.sqrt(2)*W
    M_new = Mt(a_new, eta)
    x_new = x + 0.5 * np.dot(M_new + M_old,B)

    return x_new, a_new


def create_diatom(x, l, theta):
    """ Create 3D positions from 2D center of mass, length and angle
    """
    c = np.cos(theta)
    s = np.sin(theta)
    dx  = np.array([c, s]) * l

    pos = np.zeros((2,3))
    pos[0,0:2] = x + dx
    pos[1,0:2] = x - dx
    return pos

# Parameters of the trajectory
eta     =       0  # degree of anisotropy
n_steps =  100000  # steps
angle   =       0  # initial angle in radians
com     = [0.,0.]  # initial position
length  =       1  # length of the particle in angstroms
l_box   =     100  # box dimensions in angstrom

# Create MDA Universe from scratch:
atoms    = 2
residues = 1
u = mda.Universe.empty(atoms,
                       residues,
                       atom_resindex=[0, 0],
                       trajectory=True)

u.add_TopologyAttr('names')
u.add_TopologyAttr('resids')
u.add_TopologyAttr('resnames')
u.atoms.names = ['C','C']
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
    com, angle = update_step(com, angle, eta)

# load the coordinates into the universe
# and write them in a GRO and XTC format
u.load_new(coordinates, order='fac')
u.dimensions = np.array([l_box,l_box, l_box, 90, 90, 90])
u.atoms.write('aniso.gro')
with mda.Writer("aniso.xtc", u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        u.dimensions = np.array([l_box,l_box, l_box, 90, 90, 90])
        W.write(u.atoms)
