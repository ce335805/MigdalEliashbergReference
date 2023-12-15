# Code to solve coupled gap equations

This repository contains the code to solve the coupled gap equations derived in https://arxiv.org/abs/2303.02176.

## Installation

For executing the code a `python` interpreter as well as `cmake`,a `hdf5` installation and intel's MKL are needed.
The preferred `c++` compile may be configured in the `CMakeLists.txt` file in the top directory.
Given a suitable configuration of all needed packages (finding hdf5 libraries and MKL) the installation using cmake should take less than a minute on a standard laptop.

## Running

First clone the repository and `cd` into the main directory.
Begin by creating the input data.
For this change into the `data` directory
```
$ cd data/
```
Create the directory holding the input data
```
$ mkdir input
```
and one that contains plots of inputs such as bandstructures or light-matter couplings
```
$ mkdir savedPlots
```
The bandstructure for graphene can be created by running
```
$ python createGrapheneBandstructure.py
```
This also creates the box-Coupling used as an approximation in https://arxiv.org/abs/2303.02176 which is also needed by the code.
In the repository 50x50 k-points are used which will reproduce a Tc close to the one reported in the paper which has been computed with a grid of 500x500 k-points.
The resolution may be changed in `createGrapheneBandstructure.py` by changing `N` at the beginning of the script.

For creating the input data for the sawtooth chain run

```
$ python createSawTooth.py
```
The chainlength is set to 1000 sites here again reproducing results similar to those reported in the paper where 100,000 sites were used.

Finally create the directory in which the output of the main program is saved

```
$ mkdir results
```
Change back to the top directory `$ cd ..` and install the program using `cmake`
```
$ cmake CMakeLists.txt
```
and build the project using the generated `Makefile`
```
$ make
```
In order to compute Tc for graphene for the setup proposed in https://arxiv.org/abs/2303.02176 simply run the code
```
$ ./EliashberChain.out
```

ProTipp: If the input files are not found check that line-endings are correctly configured for Linux / Windows in the input file.
In order to convert line-endings from Windows-ones to Unix-ones the tool `dos2unix` is quite convenient.

In order to reproduce the plot the computes Tc as a function of the boson frequency Omega, comment out the first part for graphene of the code and comment in the corresponding part in `src/eliashbergChain.cpp` -- see comments in that file.
Build and run the code again

```
$ make && ./EliashbergChain.out
```
Afterwards go to `parameters/globals.h` and change `drivenPhonons` to `0.`.
Build and run the code again

```
$ make && ./EliashbergChain.out
```

To produce the plot change into the plotting directory `cd plot/` and create the directory to save the plots into
```
$ mkdir savedPlots
```
and afterwards run the plotting script
```
$ python plotting.py
```
To reproduce the phase diagram some little changes are needed to the code.
First edit `parameters/globals.h` commenting out line 11 and uncommenting line 10.
Then edit `src/processing/localOneGap.cpp` uncommenting line 13 and line 456 -- also see comment in that file.
In `src/eliashbergChain.cpp` uncomment the last part for computing the phase diagram (see comment in that file) and potentially comment out the previous part on computing Tc as a function of Omega.

Build and run the code again

```
$ make && ./EliashbergChain.out
```
and produce the plot as
```
$ python plotting.py
```


