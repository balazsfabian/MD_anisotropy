# MD Anisotropy

MD Anisotropy is a collection of Python tools for
determining anistropic diffusion in Molecular Dynamics.


## Installation

In its current form, there is no need to install the package.
It consist of standalone scripts that can be executed in a workflow.


## Usage

To determine the anisotropy of a simulated particle, the following
steps must be made:

* determine rotation matrix: the rotation matrix is the operation
that must be performed on the current configuration to align it
with a reference (in the least squares sense). This can be performed
by `gmx rotmat`.

* calculate fixed (but arbitrary) angle MSD (FA-MSD): the fixed-angle
MSD calculates the full MSD Tensor relative to a fixed **initial**
angle, at lagtime = 0. There are two marginal cases of interest:
    * When the major axis of the diffusing particle in the reference
      configuration aligns with the laboratory X axis, then at small
      lagtimes the MSD along X and Y will be uncorrelated,
      with larger MSD along the X axis.
    * When the major axis is at 45 degrees, the initial displacements
      along X and Y will be completely correlated.
```bash
python diffusion_tensor.py reference.gro traj.xtc rotmat.xvg MSDtensor.dat
```

* calculate the principal axes of anisotropiy diffusion:
    The task of determining the principal axes of diffusion is then performed
    by calculating the eigenvectors and eigenvalues of the FA-MSD at a given
    lagtime. At small lagtimes, the difference between the axes is representative
    of the anisotropy, which should vanish in the long lag-time limit. However,
    because of the eigenvector determination, statistical errors are detected by
    the method and visible in the resulting curves as an apparent difference
    between the major and minor axes of diffusion.
```bash
python orient.py MSDtensor.dat MSDtensor-DIAG.dat
```
