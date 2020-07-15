# Cylindrical-Pressure-Tensor
Post-analysis code for calculating the local pressure tensor in a Lennard-Jones or molecular system having cylindrical geometries/interfaces. 

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
- 

