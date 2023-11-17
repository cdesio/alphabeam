# Lu177 Simulation

## Origin
Repo created from two repos:
  - git@github.com:lg9783/DaRT.git/savePS at commit e491462d018ba29f9ab55aab8d1c2a15dc268f7e -> decaySimulation folder
  - git@github.com:lg9783/RBE.git/inputFromDaRT at commit 420bb06e90c883b7768beacc1f064bd8fd189c93 -> DNAsimulation folder
  
## Overview
The decaySimulation simulates a point-like alpha source, shooting particles along the Z axis.
Particles are saved as they cross a voxelised geometry.
The produced output files are:
- an output.root file, containing basic info, e.g. number of primaries
- an output.bin file, to be used as input for the DNA simulation (RBE), with the option -PS

## How to Run

./alphaBeam -mac alphaBeam.in -out output (optional: -gui)


## Other info

The default macro generates 5.5 MeV alphas along the Z axis.
To change the number of alphas, change the value to /run/beamOn 

Voxels are placed in src/DetectorConstruction.cc
Default number of voxels is 12000: 20 along x, 10 along y and 60 along z.

A "copy number" is assigned to each voxel, identifying the layer across the Z axis (all voxels with the same Z have the same copyNo.). This is used in the RBE clustering to group the events and evaluate the DNA damage as a function of the distance. 
