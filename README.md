# DimerGen

DimerGen is a molecular utility program that can be used to generate dimers in a 6D space. The dimers are
described by the center-of-mass (COM) distance between their constituent monomers and 5 Euler angles 
(in the ZYZ convention). The user must provide the geometry of each monomer (in XYZ file format), along 
with the 6 coordinates: $( R_{\rm COM}, \beta_A, \gamma_A,  \alpha_B, \beta_B, \gamma_B )$. Note: $\alpha_A$ is hardcoded and set to 0, so it is not a coordinate.

There are two versions of the program, one written in C++ and another written in Python. The former can be
found in the [cppDimerGen](https://github.com/jwmelkumov/DimerGen/tree/main/cppDimerGen) subdirectory and the latter can be found in the [pyDimerGen](https://github.com/jwmelkumov/DimerGen/tree/main/pyDimerGen) subdirectory.

## Table of Contents  
- [Euler Angles and Rigid Body Rotations](#euler-angles-and-rigid-body-rotations) 
- [Notes](#notes) 
- [Dependencies](#dependencies) 
- [Usage](#usage) 

## Euler Angles and Rigid Body Rotations

<p align="center">
<img src="https://raw.githubusercontent.com/jwmelkumov/DimerGen/main/Documentation/EulerRotations.png" width="500">
</p>

Given a set of atomic coordinates and masses (automatically extracted via lookup table from the substrings of
atomic labels in the input XYZ file) for each monomer, the center-of-mass is computed and the coordinates are shifted such that the COM is at the origin. The inertia tensor is then computed for each monomer and the principal axes are obtained by solving for the eigenvectors of the inertia tensor. After the principal axes are obtained, the coordinates of each monomer are transformed into the principal axes frame. Using the user-specified Euler angles, the rotation matrix is constructed and applied to the coordinates in the principal axes frame. The rotated coordinates are then transformed back to the original frame and the separation vector between the two monomers is used to translate the coordinates by the user-specified separation (in angstroms).

The convention used here is ZYZ, so the rotation matrix is constructed as:
<p align="center">
<img src="https://raw.githubusercontent.com/jwmelkumov/DimerGen/main/Documentation/R.png" width="500">
</p>

where
<p align="center">
<img src="https://raw.githubusercontent.com/jwmelkumov/DimerGen/main/Documentation/Rz.png" width="500">
</p>

and
<p align="center">
<img src="https://raw.githubusercontent.com/jwmelkumov/DimerGen/main/Documentation/Ry.png" width="500">
</p>

## Notes
- DimerGen supports monomers with off-atomic sites (such as M in TIP4P and LP "lone pairs" in force fields).
- Monomers containing all elements are supported.
- This is not a simulation package. As such, there is no potential energy function or scoring function implemented. This is a utility that is meant to be paired with external automation and logic.
## Dependencies
- Eigen3 (for the C++ version)
- C++ compiler (e.g., g++, for the C++ version)
- NumPy (for the Python version)
- Python3 (for the Python version)
## Usage 
To compile the C++ program, one can use the following line:
```bash
$ g++ -o DimerGen -I <path_to_eigen3> DimerGen.cpp
```
Once the code has been compiled, one can then run the program with the following line:
```bash
./DimerGen xyzfileA xyzfileB sep betaA gammaA alphaB betaB gammaB
```
With the Python program, one can achieve the same result by running the program with the following line:
```bash
./pyDimerGen.py xyzfileA xyzfileB sep betaA gammaA alphaB betaB gammaB
```
