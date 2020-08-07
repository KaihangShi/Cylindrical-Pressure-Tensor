# Cylindrical-Pressure-Tensor
Post-analysis code for the local pressure tensor in a Lennard-Jones or molecular system having cylindrical geometries/interfaces. 

Reference: Shi et al. https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00607

## Requirement
gfortran

## Files
For Lennard-Jones system:

/Lennard_Jones_fluid/
- cylinpress_IK_ljcut.f90 (local pressure tensor using Irving-Kirkwood definition)
- cylinpress_H_ljcut.f90 (local pressure tensor using Harasima definition)

For molecular system (demonstration of rigid water in single-wall carbon nanotube):

/Water_SWCNT/
- cylinpress_IK_ljcoulcut.f90 (local axial pressure using Irving-Kirkwood definition with simple cutoff)

/Water_SWCNT/Harasima_Ewald/
- cylinpress_H_ewald.f90 (main program; local axial pressure using Harasima/Ewald method)
- ewald_mol.f90 (subroutine to evaluate pressure in reciprocal space)
- global.f90 (module containing global data)
- set_ewld.f90 (subroutine to set up Ewald parameters)

## Usage
1. Set up the system parameters (temperature, potential paramters, profile resolution etc.) at the beginning of the main program.
2. compile the fortran code as usual.
