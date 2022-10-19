# Readdata.f
## Authors: Ashley Liang, Hans Mueller

## Overview
This fortran file was part of the computational astronomy research project I conducted under Professor Hans Mueller. This file's purpose was to take the data file produced by the existing ZEUS-3D program, take the necessary data from it to perform calculations to factor frequency in (as radial velocity is very important in astronomy) and get a new value dtau (the overall amount of light that passes through at certain frequencies/wavelengths) which will be printed out to a dtau data file and absorb data file. More information about each of the files can be found in readdata.f. 

## My experience
Creating this program was very challenging for me, both because I've never learned Fortran before and because I had a hard time grasping the physics behind it and applying it into the code. My professor started the program and coded the first two loops that reads through the files, stores them into arrays, and derives the first few variables. I developed everything after that. My professor sent me readings with equations that I needed to use in the program, but it was pretty difficult for me to understand which data points from the data file corresponded to which variable in the equations.

## Usage
To compile, run `gfortran -c readdata.f wofz.f` and then `gfortran *.o`. To execute the code, run:
```
./a.out filename dtaufile absorbfile
```
For example:
```
./a.out lastframem5v100.dat dtau5v100.dat absorb5v100.dat
```

## Results from Data


